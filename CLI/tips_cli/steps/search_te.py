from __future__ import annotations
import shlex
import ast
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional
import pandas as pd
import yaml,os
from tips_cli.runner import DockerRunner


@dataclass
class SearchTeOutputs:
    sample_path: Path
    search_dir: Path
    mzml_dir: Path
    log_file: Path
    commands_sh: Path


# ----------------------------
# Basic FS helpers
# ----------------------------
def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _nonempty(p: Path) -> bool:
    return p.exists() and p.is_file() and p.stat().st_size > 0


def _write_cmd(commands_sh: Path, cmd_line: str) -> None:
    _ensure_dir(commands_sh.parent)
    with commands_sh.open("a", encoding="utf-8") as f:
        f.write(cmd_line.rstrip() + "\n")


def _safe_symlink(src: Path, dst: Path) -> None:
    """Create a symlink if dst does not exist; overwrite only if dst is a broken link."""
    _ensure_dir(dst.parent)
    if dst.exists():
        return
    if dst.is_symlink() and not dst.resolve().exists():
        dst.unlink()
    dst.symlink_to(src)


def _collect_extra_mounts(mount_root: Path, paths: List[Path]) -> List[Path]:
    """Collect parent directories that must be bind-mounted for docker to resolve symlink targets.

    If a resource file is outside mount_root, the bind mount of mount_root alone is insufficient.
    We add the resource's parent directory as an additional same-path mount.
    """
    extra: List[Path] = []
    seen: set[str] = set()
    for p in paths:
        p = p.resolve()
        try:
            p.relative_to(mount_root)
            continue
        except Exception:
            parent = p.parent
            key = str(parent)
            if key not in seen:
                seen.add(key)
                extra.append(parent)
    return extra

def _collect_blastdb_target_mounts(
    mount_root: Path,
    blastdb_prefix_in_sample: Path,
) -> List[Path]:
    """Collect extra bind-mount dirs for BLAST DB symlink targets.

    We expect prefix.{pin,psq,phr} to exist under the sample tree, but they may be symlinks
    pointing outside mount_root (e.g. /data/.../genome/...).
    BLAST resolves symlinks and reads the real files, so those real parent dirs must be mounted.
    """
    extra: List[Path] = []
    seen: set[str] = set()

    for suf in (".pin", ".psq", ".phr"):
        p = blastdb_prefix_in_sample.with_suffix(suf)
        if not p.exists():
            continue
        if not p.is_symlink():
            continue

        target = p.resolve()
        parent = target.parent
        try:
            parent.relative_to(mount_root.resolve())
            continue
        except Exception:
            key = str(parent)
            if key not in seen:
                seen.add(key)
                extra.append(parent)

    return extra


def _load_cfg(cfg_path: Path) -> Dict:
    cfg = yaml.safe_load(cfg_path.read_text(encoding="utf-8")) or {}
    return cfg


# ----------------------------
# FASTA helpers (build combined target/decoy db)
# ----------------------------
def _iter_fasta_records(fa_path: Path) -> Iterable[tuple[str, str]]:
    """Yield (header, sequence) for a FASTA file (single- or multi-line sequences)."""
    header: Optional[str] = None
    seq_chunks: List[str] = []
    with fa_path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            yield header, "".join(seq_chunks)


def _write_fasta(records: Iterable[tuple[str, str]], out_path: Path) -> None:
    _ensure_dir(out_path.parent)
    with out_path.open("w", encoding="utf-8") as out:
        for h, s in records:
            out.write(f">{h}\n{s}\n")


def _build_combined_db(
    te_fasta: Path,
    human_fasta: Path,
    contaminants_fasta: Path,
    out_target: Path,
    out_target_decoy: Path,
) -> None:
    """Merge TE + canonical + contaminants and append reversed decoys with 'rev_' prefix."""

    def target_iter() -> Iterable[tuple[str, str]]:
        yield from _iter_fasta_records(human_fasta)
        # Keep the same concatenation order as the original TIPs pipeline:
        # canonical (human) -> TE (denovo-derived) -> contaminants
        yield from _iter_fasta_records(te_fasta)
        yield from _iter_fasta_records(contaminants_fasta)

    _write_fasta(target_iter(), out_target)

    def target_decoy_iter() -> Iterable[tuple[str, str]]:
        for h, s in _iter_fasta_records(out_target):
            yield h, s
        for h, s in _iter_fasta_records(out_target):
            yield f"rev_{h}", s[::-1]

    _write_fasta(target_decoy_iter(), out_target_decoy)


# ----------------------------
# Params helpers (MSFragger runtime params)
# ----------------------------
def _rewrite_kv_params_file(
    template_path: Path,
    out_path: Path,
    overrides: Dict[str, str],
) -> Path:
    """Rewrite a simple key=value params file with overrides.

    Used for MSFragger params. Updates existing keys if found, otherwise appends missing keys.
    """
    _ensure_dir(out_path.parent)

    lines = template_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    key_found = {k: False for k in overrides.keys()}

    new_lines: List[str] = []
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            new_lines.append(line)
            continue

        if "=" in stripped:
            k = stripped.split("=", 1)[0].strip()
            if k in overrides:
                new_lines.append(f"{k} = {overrides[k]}")
                key_found[k] = True
                continue

        new_lines.append(line)

    for k, v in overrides.items():
        if not key_found[k]:
            new_lines.append(f"{k} = {v}")

    out_path.write_text("\n".join(new_lines) + "\n", encoding="utf-8")
    return out_path


# ----------------------------
# pepXML helpers
# ----------------------------
def _rename_pepxml_to_pepxml_xml(work_dir: Path) -> None:
    """Normalize pepXML extensions to .pep.xml (needed for downstream TPP scripts)."""
    for p in work_dir.iterdir():
        if p.is_file() and p.name.endswith(".pepXML"):
            p.rename(p.with_name(p.name.replace(".pepXML", ".pep.xml")))



def _collect_engine_outputs_for_mzml(mzml: Path, raw_dir: Path) -> None:
    """Move all result files produced for a given mzML into raw_dir.

    Comet often writes outputs next to the input mzML. We collect all files whose names
    start with mzml.stem from mzML directory and move them into raw_dir.
    """
    src_dir = mzml.parent
    prefix = mzml.stem

    candidates = [p for p in src_dir.glob(f"{prefix}*") if p.is_file()]
    if not candidates:
        return

    _ensure_dir(raw_dir)
    for p in candidates:
        if p.resolve() == mzml.resolve():
            continue

        dst = raw_dir / p.name
        if dst.exists():
            dst = raw_dir / f"{p.stem}__dup{p.suffix}"

        shutil.move(str(p), str(dst))

    _rename_pepxml_to_pepxml_xml(raw_dir)


def _split_pepxml_te_only_in_container(
    runner: DockerRunner,
    commands_sh: Path,
    log_file: Path,
    te_fasta: Path,
) -> None:
    """Split pepXML into TE-only pepXML files inside container (python+lxml).

    This mirrors the original pipeline behavior: keep only search_hit where
    'Denovo' appears in search_hit@protein, and drop spectrum_query if no hits remain.
    """
    script = r"""
import glob
import os
from lxml import etree

TE_FASTA = os.environ.get("TE_FASTA", "").strip()

def load_te_headers(te_fasta_path: str) -> set[str]:
    headers = set()
    if not te_fasta_path or not os.path.exists(te_fasta_path):
        return headers
    with open(te_fasta_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                h = line[1:].strip()
                if not h:
                    continue
                headers.add(h)
                headers.add(h.split()[0])
    return headers

TE_HEADERS = load_te_headers(TE_FASTA)


def normalize_ext():
    for p in glob.glob("*.pepXML"):
        os.rename(p, p.replace(".pepXML", ".pep.xml"))

def is_te_hit(protein: str) -> bool:
    if not protein:
        return False

    # Protein strings in pepXML may contain extra annotations after '#' or whitespace.
    raw = str(protein)
    p0 = raw.split("#")[0].strip()
    p1 = p0.split()[0] if p0 else ""

    # First try exact matching against TE FASTA headers if available.
    if TE_HEADERS and ((p0 in TE_HEADERS) or (p1 in TE_HEADERS)):
        return True

    # Fallback rule aligned with the original TIPs pipeline:
    # keep hits only if the pepXML protein string contains the token 'Denovo'.
    # This avoids incorrectly classifying canonical proteins as TE due to substring matches.
    return ("Denovo" in raw)


def split_one(in_path: str) -> str:
    parser = etree.XMLParser(remove_blank_text=False, recover=True, huge_tree=True)
    tree = etree.parse(in_path, parser)
    root = tree.getroot()

    # root -> msms_run_summary -> spectrum_query -> search_result -> search_hit
    # Namespace-safe checks: we use 'endswith' on tag localname.
    for run_summary in list(root):
        if not run_summary.tag.endswith("msms_run_summary"):
            continue

        for spectrum_query in list(run_summary):
            if not spectrum_query.tag.endswith("spectrum_query"):
                continue
            if len(spectrum_query) == 0:
                continue

            search_result = spectrum_query[0]
            # Remove non-TE hits
            for hit in list(search_result):
                if not hit.tag.endswith("search_hit"):
                    continue
                protein = hit.get("protein", "")
                if not is_te_hit(protein):
                    search_result.remove(hit)

            # Drop spectrum_query if no hits remain
            #if len(search_result) == 0:
                # run_summary.remove(spectrum_query)

    out_path = in_path.replace(".pep.xml", "_TE.pep.xml")
    tree.write(out_path, xml_declaration=True, encoding="UTF-8", pretty_print=True)
    return out_path

def main():
    normalize_ext()
    files = sorted(glob.glob("*.pep.xml"))
    if not files:
        raise SystemExit("No .pep.xml under workdir")
    for f in files:
        if f.endswith("_TE.pep.xml"):
            continue
        split_one(f)

if __name__ == "__main__":
    main()
"""
    inner_bash = [
        "bash",
        "-lc",
        f"export TE_FASTA={shlex.quote(str(te_fasta))}\npython - << 'PY'\n{script}\nPY",
    ]
    cmd_line = runner.cmd_to_shell(runner.build_run_cmd(inner_bash))
    _write_cmd(commands_sh, cmd_line)
    rc = runner.run(inner_bash, dry_run=False, log_file=log_file)
    if rc != 0:
        raise RuntimeError(f"Split pepXML TE-only failed (exit {rc})")


# ----------------------------
# FDR and blastp removal (TE-only)
# ----------------------------
def _filter_by_fdr(pepxml_csv: Path, fdr_threshold: float) -> Path:
    """Replicate original pipeline filter_by_fdr (PSM-level, charge-stratified).

    This function must be robust to empty/invalid pepxml2csv outputs (0 hits),
    in which case it writes an empty TSV with a header and returns.
    """
    out_path = pepxml_csv.with_name(pepxml_csv.name.replace(".pep.csv", "_fdr.csv"))

    # 1) Read input safely
    try:
        df_all = pd.read_csv(pepxml_csv, sep="\t")
    except pd.errors.EmptyDataError:
        pd.DataFrame({"peptide": []}).to_csv(out_path, index=False, sep="\t")
        return out_path
    except Exception:
        pd.DataFrame({"peptide": []}).to_csv(out_path, index=False, sep="\t")
        return out_path

    # If the file has no columns or no rows, return an empty table with a reasonable header.
    if df_all is None or df_all.empty or df_all.shape[1] == 0:
        pd.DataFrame({"peptide": []}).to_csv(out_path, index=False, sep="\t")
        return out_path

    # Keep an empty template (same columns) for consistent downstream parsing.
    empty_template = df_all.head(0).copy()

    # 2) Validate required columns; if missing, return empty template
    required_cols = {"protein", "analysis_result", "assumed_charge"}
    if not required_cols.issubset(set(df_all.columns)):
        empty_template.to_csv(out_path, index=False, sep="\t")
        return out_path

    # 3) If too few PSMs, skip FDR model (heuristic guard for small tables)
    if len(df_all) < 20:
        empty_template.to_csv(out_path, index=False, sep="\t")
        return out_path

    # 4) Compute peptideprophet probability and decoy flag
    df_all["protein_primary"] = df_all["protein"].apply(lambda x: str(x).split("#")[0])

    try:
        df_all["peptideProphet"] = df_all["analysis_result"].apply(ast.literal_eval)
        df_all["peptideProphet_probability"] = df_all["peptideProphet"].apply(
            lambda x: float(x[0]["peptideprophet_result"]["probability"])
        )
    except Exception:
        # If analysis_result cannot be parsed, return empty.
        empty_template.to_csv(out_path, index=False, sep="\t")
        return out_path

    df_all["decoy"] = df_all["protein_primary"].apply(lambda x: 1 if "rev_" in str(x) else 0)
    df_all = df_all.sort_values(by=["peptideProphet_probability"], ascending=False).reset_index(drop=True)

    # 5) Charge-stratified FDR
    df_merged = pd.DataFrame()
    for charge, _ in df_all.groupby("assumed_charge"):
        df = df_all[df_all["assumed_charge"] == charge].copy()

        # Avoid division by zero in pathological cases
        df["target_cumsum"] = (df["decoy"] == 0).cumsum()
        df["decoy_cumsum"] = (df["decoy"] == 1).cumsum()
        df = df[df["target_cumsum"] > 0].copy()
        if df.empty:
            continue

        df["fdr"] = df["decoy_cumsum"] / df["target_cumsum"]
        df.reset_index(drop=True, inplace=True)

        reversed_df = df.iloc[::-1]
        ok = reversed_df[reversed_df["fdr"] < fdr_threshold]
        if len(ok) == 0:
            continue

        first_index = ok.index[0]
        filtered_df = df.iloc[: first_index + 1]
        filtered_df = filtered_df[filtered_df["decoy"] == 0]

        if not filtered_df.empty:
            df_merged = filtered_df if df_merged.empty else pd.concat([df_merged, filtered_df], ignore_index=True)

    # 6) Write output
    if df_merged.empty:
        empty_template.to_csv(out_path, index=False, sep="\t")
        return out_path

    df_merged = df_merged.sort_values(by=["peptideProphet_probability"], ascending=False)
    df_merged.to_csv(out_path, index=False, sep="\t")
    return out_path

def _blastp_remove_canonical(
    runner: DockerRunner,
    commands_sh: Path,
    log_file: Path,
    te_fdr_csv: Path,
    human_blastdb_prefix: Path,
    blastp_bin: str,
    task: str,
    evalue: float,
    outfmt: str,
    threads: int,
) -> Path:
    """Run blastp-short against canonical proteome and remove full-length 100% identity matches."""
    # Fast-path: skip BLASTP if input is empty or missing.
    if (not te_fdr_csv.exists()) or te_fdr_csv.stat().st_size == 0:
        out_path = te_fdr_csv.with_suffix("").with_name(te_fdr_csv.stem + "_blastp.filter.txt")
        pd.DataFrame().to_csv(out_path, index=False, sep="\t")
        return out_path

    try:
        df = pd.read_csv(te_fdr_csv, sep="\t")
    except pd.errors.EmptyDataError:
        out_path = te_fdr_csv.with_suffix("").with_name(te_fdr_csv.stem + "_blastp.filter.txt")
        pd.DataFrame().to_csv(out_path, index=False, sep="\t")
        return out_path

    # If the table has no rows, nothing to BLAST.
    if df.empty:
        out_path = te_fdr_csv.with_suffix("").with_name(te_fdr_csv.stem + "_blastp.filter.txt")
        df.to_csv(out_path, index=False, sep="\t")
        return out_path

    peptide_col = "peptide" if "peptide" in df.columns else ("Peptide" if "Peptide" in df.columns else None)

    if peptide_col is None:
        raise ValueError("TE FDR CSV must contain 'peptide' or 'Peptide' column")

    query_fasta = te_fdr_csv.with_suffix("").with_name(te_fdr_csv.stem + "_blastp.fasta")
    with query_fasta.open("w", encoding="utf-8") as f:
        for seq in sorted(set(df[peptide_col].astype(str).tolist())):
            f.write(f">{seq}\n{seq}\n")

    blast_out = te_fdr_csv.with_suffix("").with_name(te_fdr_csv.stem + "_blastpResult.txt")
    inner = [
        blastp_bin,
        "-task",
        task,
        "-query",
        str(query_fasta),
        "-db",
        str(human_blastdb_prefix),
        "-out",
        str(blast_out),
        "-outfmt",
        outfmt,
        "-evalue",
        str(evalue),
        "-num_threads",
        str(threads),
    ]
    cmd_line = runner.cmd_to_shell(runner.build_run_cmd(inner))
    _write_cmd(commands_sh, cmd_line)
    rc = runner.run(inner, dry_run=False, log_file=log_file)
    if rc != 0:
        raise RuntimeError(f"blastp failed with exit code {rc}")

    cols = outfmt.split()[1:]
    bdf = pd.read_csv(blast_out, sep="\t", header=None, names=cols)
    if "qaccver" not in bdf.columns or "pident" not in bdf.columns:
        raise ValueError("blastp outfmt must include qaccver and pident")

    bdf["seq_len"] = bdf["qaccver"].astype(str).apply(len)
    if "qstart" in bdf.columns and "qend" in bdf.columns:
        bdf = bdf[(bdf["qstart"] == 1) & (bdf["qend"] == bdf["seq_len"])]
    bdf = bdf[bdf["pident"] == 100]
    remove_peptides = set(bdf["qaccver"].astype(str).unique().tolist())

    out_df = df[~df[peptide_col].astype(str).isin(remove_peptides)]
    out_path = te_fdr_csv.with_suffix("").with_name(te_fdr_csv.stem + "_blastp.filter.txt")
    out_df.to_csv(out_path, index=False, sep="\t")
    return out_path

def _blastdb_exists(prefix: Path) -> bool:
    """Return True if BLAST protein DB index files exist for given prefix."""
    return (
        prefix.with_suffix(".pin").exists()
        and prefix.with_suffix(".psq").exists()
        and prefix.with_suffix(".phr").exists()
    )

def _ensure_human_blastdb(
    runner: DockerRunner,
    commands_sh: Path,
    log_file: Path,
    human_fasta: Path,
    requested_prefix: Optional[str],
    out_dir: Path,
) -> Path:
    """Ensure a canonical proteome BLAST DB exists; build it with makeblastdb if needed.

    - If requested_prefix is relative, resolve it under runner.mount_root (sample root).
    - If index files are missing, run makeblastdb to create them.
    - Always return an absolute path usable inside the container.
    """
    # Decide prefix path
    if requested_prefix and str(requested_prefix).strip():
        prefix = Path(str(requested_prefix).strip())
        # If user config gives a relative path, treat it as relative to the sample root
        if not prefix.is_absolute():
            prefix = runner.mount_root / prefix
    else:
        prefix = out_dir / "human_proteome"

    prefix = prefix.resolve()
    _ensure_dir(prefix.parent)

    # Build DB if missing
    if not _blastdb_exists(prefix):
        inner = [
            "makeblastdb",
            "-in",
            str(human_fasta),
            "-dbtype",
            "prot",
            "-out",
            str(prefix),
        ]
        cmd_line = runner.cmd_to_shell(runner.build_run_cmd(inner))
        _write_cmd(commands_sh, cmd_line)
        rc = runner.run(inner, dry_run=False, log_file=log_file)
        if rc != 0:
            raise RuntimeError(f"makeblastdb failed with exit code {rc} (prefix={prefix})")

    return prefix



# ----------------------------
# Main step
# ----------------------------
def search_te_step(
    outputs: SearchTeOutputs,
    image: str,
    config_path: Path,
    engines: List[str],
    mzml_glob: str,
    fdr: float,
    blastp_threads: int,
    dry_run: bool,
    resume: bool,
    overwrite: bool,
) -> None:
    """Run one or more DB search engines and perform TE-only xinteract + PSM FDR + blastp removal.

    Notes:
    - No multi-engine merging in this step.
    - TE-only: we split pepXML into *_TE.pep.xml inside container (python+lxml),
      then run xinteract only on TE pepXMLs.
    """
    _ensure_dir(outputs.search_dir)
    _ensure_dir(outputs.log_file.parent)
    _ensure_dir(outputs.commands_sh.parent)

    cfg = _load_cfg(config_path)

    # Resources: symlink user-provided fasta/params into sample/resources/
    resources_dir = outputs.sample_path / "resources"
    res_db_dir = resources_dir / "db"
    res_params_dir = resources_dir / "search_params"
    _ensure_dir(res_db_dir)
    _ensure_dir(res_params_dir)

    human_fa = Path(str(cfg.get("db", {}).get("human_fasta", "") or "")).expanduser()
    contaminants_fa = Path(str(cfg.get("db", {}).get("contaminants_fasta", "") or "")).expanduser()
    if not str(human_fa).strip():
        raise ValueError("db.human_fasta is empty in config")
    if not str(contaminants_fa).strip():
        raise ValueError("db.contaminants_fasta is empty in config")

    # Resolve relative paths against the sample root (not the current working directory)
    if not human_fa.is_absolute():
        human_fa = outputs.sample_path / human_fa
    if not contaminants_fa.is_absolute():
        contaminants_fa = outputs.sample_path / contaminants_fa

    human_fa = human_fa.resolve()
    contaminants_fa = contaminants_fa.resolve()
    if not human_fa.exists():
        raise FileNotFoundError(f"human fasta not found: {human_fa}")
    if not contaminants_fa.exists():
        raise FileNotFoundError(f"contaminants fasta not found: {contaminants_fa}")


    human_fa = human_fa.resolve()
    contaminants_fa = contaminants_fa.resolve()
    if not human_fa.exists():
        raise FileNotFoundError(f"human fasta not found: {human_fa}")
    if not contaminants_fa.exists():
        raise FileNotFoundError(f"contaminants fasta not found: {contaminants_fa}")

    human_link = res_db_dir / human_fa.name
    cont_link = res_db_dir / contaminants_fa.name
    if overwrite:
        if human_link.is_symlink() or human_link.exists():
            human_link.unlink(missing_ok=True)
        if cont_link.is_symlink() or cont_link.exists():
            cont_link.unlink(missing_ok=True)
    _safe_symlink(human_fa, human_link)
    _safe_symlink(contaminants_fa, cont_link)

    params_cfg = cfg.get("search_params", {}) or {}
    param_src_paths: List[Path] = []
    for k in ("comet", "msfragger", "msgfplus"):
        v = str(params_cfg.get(k, "") or "").strip()
        if not v:
            continue
        src = Path(v).expanduser()
        # Resolve relative paths against the sample root (not the current working directory)
        if not src.is_absolute():
            src = outputs.sample_path / src
        src = src.resolve()

        if not src.exists():
            raise FileNotFoundError(f"search_params.{k} not found: {src}")
        param_src_paths.append(src)
        dst = res_params_dir / src.name
        if overwrite and (dst.is_symlink() or dst.exists()):
            dst.unlink(missing_ok=True)
        _safe_symlink(src, dst)

    # Combined DB
    te_fasta = outputs.sample_path / "denovo" / "te" / "Denovo_TE_SoftMerged_Merged.fasta"
    if not te_fasta.exists():
        raise FileNotFoundError(f"TE fasta not found: {te_fasta} (run denovo-te first)")

    db_dir = outputs.search_dir / "db"
    _ensure_dir(db_dir)
    combined_target = db_dir / "combined.fasta"
    combined_target_decoy = db_dir / "combined_target_decoy.fasta"
    if overwrite or not (combined_target.exists() and combined_target_decoy.exists()):
        _build_combined_db(te_fasta, human_link, cont_link, combined_target, combined_target_decoy)

    # Inputs
    mzmls = sorted(outputs.mzml_dir.glob(mzml_glob))
    if not mzmls:
        raise FileNotFoundError(f"No mzML found: {outputs.mzml_dir}/{mzml_glob}")

    # Base mounts needed for resolving symlinked FASTA/params outside the sample root
    extra_mounts = _collect_extra_mounts(outputs.sample_path, [human_fa, contaminants_fa] + param_src_paths)

    # If user provides a BLAST DB prefix (B plan), it is typically under sample/resources,
    # but the actual .pin/.psq/.phr may be symlinks pointing outside and must be mounted too.
    human_blast_prefix_cfg = str(cfg.get("db", {}).get("human_blastdb_prefix", "") or "").strip()
    if human_blast_prefix_cfg:
        prefix_in_sample = Path(human_blast_prefix_cfg)
        if not prefix_in_sample.is_absolute():
            prefix_in_sample = outputs.sample_path / prefix_in_sample
        # Add mounts for symlink targets of prefix.{pin,psq,phr}
        extra_mounts.extend(_collect_blastdb_target_mounts(outputs.sample_path, prefix_in_sample))

    # De-duplicate mounts (keep order)
    dedup: List[Path] = []
    seen: set[str] = set()
    for p in extra_mounts:
        key = str(p)
        if key not in seen:
            seen.add(key)
            dedup.append(p)
    extra_mounts = dedup

    runner = DockerRunner(image=image, mount_root=outputs.sample_path, user_map=True, extra_mounts=extra_mounts)

    # Canonical blast db (for TE-only removal)
    blastdb_out_dir = res_db_dir / "blastdb"
    _ensure_dir(blastdb_out_dir)
    human_blastdb_prefix = _ensure_human_blastdb(
        runner=runner,
        commands_sh=outputs.commands_sh,
        log_file=outputs.log_file,
        human_fasta=human_link,
        requested_prefix=human_blast_prefix_cfg if human_blast_prefix_cfg else None,
        out_dir=blastdb_out_dir,
    )

    blast_cfg = cfg.get("blastp", {}) or {}
    blastp_bin = str(blast_cfg.get("bin", "blastp"))
    blastp_task = str(blast_cfg.get("task", "blastp-short"))
    blastp_evalue = float(blast_cfg.get("evalue", 20000))
    blastp_outfmt = str(
        blast_cfg.get(
            "outfmt",
            "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore",
        )
    )

    # Per-engine execution
    for engine in engines:
        if engine not in {"comet", "msfragger", "msgfplus"}:
            raise ValueError(f"Unsupported engine: {engine}")

        engine_dir = outputs.search_dir / engine
        raw_dir = engine_dir / "raw"
        te_dir = engine_dir / "te"

        final_marker = te_dir / f"{engine}_TE_merged_fdr_blastp.filter.txt"
        if resume and _nonempty(final_marker) and not overwrite:
            continue

        if overwrite and engine_dir.exists():
            shutil.rmtree(engine_dir)

        _ensure_dir(raw_dir)
        _ensure_dir(te_dir)

        engine_runner = DockerRunner(
            image=image,
            mount_root=outputs.sample_path,
            user_map=True,
            workdir_in_container=str(raw_dir),
            extra_mounts=extra_mounts,
        )

        # 1) DB search
        if engine == "comet":
            comet_params = str(params_cfg.get("comet", "") or "").strip()
            if not comet_params:
                raise ValueError("search_params.comet is empty in config")

            comet_param_path = (res_params_dir / Path(comet_params).name).resolve()

            link_paths = []
            for mzml in mzmls:
                link_path = raw_dir / mzml.name
                if not link_path.exists() and not link_path.is_symlink():
                    os.symlink(mzml, link_path)
                link_paths.append(str(link_path))

            inner = ["/opt/tpp/bin/comet", f"-P{comet_param_path}", f"-D{combined_target_decoy}"] + link_paths
            cmd_line = engine_runner.cmd_to_shell(engine_runner.build_run_cmd(inner))
            _write_cmd(outputs.commands_sh, cmd_line)

            if not dry_run:
                rc = engine_runner.run(inner, dry_run=False, log_file=outputs.log_file)
                if rc != 0:
                    raise RuntimeError(f"Comet failed (exit {rc})")

                _rename_pepxml_to_pepxml_xml(raw_dir)
            # =================================================


        elif engine == "msfragger":

            fr_params = str(params_cfg.get("msfragger", "") or "").strip()

            if not fr_params:
                raise ValueError("search_params.msfragger is empty in config")

            fr_template_path = (res_params_dir / Path(fr_params).name).resolve()

            fr_runtime_path = raw_dir / "msfragger.runtime.params"


            _rewrite_kv_params_file(

                template_path=fr_template_path,

                out_path=fr_runtime_path,

                overrides={

                    "database_name": str(combined_target_decoy),

                },

            )

            # ================= 核心修改区域 =================

            # 2. 批量建立软链接，并将所有新路径存入列表

            link_paths = []

            for mzml in mzmls:

                link_path = raw_dir / mzml.name

                if not link_path.exists() and not link_path.is_symlink():
                    os.symlink(mzml, link_path)

                link_paths.append(str(link_path))


            inner = ["java", "-Xmx16G", "-jar", "/opt/msfragger/MSFragger-4.1.jar", str(fr_runtime_path)] + link_paths

            cmd_line = engine_runner.cmd_to_shell(engine_runner.build_run_cmd(inner))

            _write_cmd(outputs.commands_sh, cmd_line)

            if not dry_run:

                rc = engine_runner.run(inner, dry_run=False, log_file=outputs.log_file)

                if rc != 0:
                    raise RuntimeError(f"MSFragger failed (exit {rc})")

                _rename_pepxml_to_pepxml_xml(raw_dir)

            # =================================================


        elif engine == "msgfplus":
            for mzml in mzmls:
                link_path = raw_dir / mzml.name
                if link_path.exists() or link_path.is_symlink():
                    continue
                os.symlink(mzml, link_path)

            mg_conf = str(params_cfg.get("msgfplus", "") or "").strip()

            if not mg_conf:
                raise ValueError("search_params.msgfplus is empty in config")

            mg_conf_path = (res_params_dir / Path(mg_conf).name).resolve()

            # MSGF+ builds a .cnlcp index for the FASTA database. When multiple processes start
            # simultaneously and the index does not exist yet, index creation can race and cause
            # intermittent failures or missing outputs. The original TIPs pipeline mitigates this
            # by starting the first task and giving it a head start before launching the rest.

            q = shlex.quote

            def build_one_job(mzml_path: Path) -> str:
                mzid_out = raw_dir / f"{mzml_path.stem}.mzid"
                pepxml_out = raw_dir / f"{mzml_path.stem}.pepXML"
                cmd1 = (
                    f"java -jar /opt/msgf/MSGFPlus.jar "
                    f"-conf {q(str(mg_conf_path))} "
                    f"-d {q(str(combined_target_decoy))} "
                    f"-s {q(str(mzml_path))} "
                    f"-o {q(str(mzid_out))}"
                )
                # mzid_conver.sh produces a modified mzid whose name is derived from pepxml_out
                mod_mzid = raw_dir / f"{mzml_path.stem}.pepXML_modify.mzid"
                cmd2 = f"bash /opt/scripts/mzid_conver.sh {q(str(mzid_out))} {q(str(pepxml_out))}"
                cmd3 = f"/opt/pwiz/idconvert {q(str(mod_mzid))} --pepXML -o {q(str(raw_dir))}"
                return f"{cmd1} && {cmd2} && {cmd3}"

            cnlcp_path = Path(str(combined_target_decoy)).with_suffix(".cnlcp")

            # If index does not exist, run the first mzML serially to create it, then run the rest in parallel.
            remaining = list(mzmls)
            if (not dry_run) and remaining and (not cnlcp_path.exists()):
                first = remaining.pop(0)
                first_bash = ["bash", "-lc", "set -euo pipefail\n" + build_one_job(first)]
                first_cmd_line = engine_runner.cmd_to_shell(engine_runner.build_run_cmd(first_bash))
                _write_cmd(outputs.commands_sh, first_cmd_line)
                rc = engine_runner.run(first_bash, dry_run=False, log_file=outputs.log_file)
                if rc != 0:
                    raise RuntimeError(f"MSGF+ failed for {first} (exit {rc})")

            # Run remaining mzMLs in parallel.
            if remaining:
                bash_lines: List[str] = ["set -euo pipefail"]
                for mzml in remaining:
                    bash_lines.append(f"( {build_one_job(mzml)} ) &")
                bash_lines.append("wait")
                inner_bash = ["bash", "-lc", "\n".join(bash_lines)]
                cmd_line = engine_runner.cmd_to_shell(engine_runner.build_run_cmd(inner_bash))
                _write_cmd(outputs.commands_sh, cmd_line)
                if not dry_run:
                    rc = engine_runner.run(inner_bash, dry_run=False, log_file=outputs.log_file)
                    if rc != 0:
                        raise RuntimeError(f"MSGF+ parallel run failed (exit {rc})")

            if not dry_run:
                _rename_pepxml_to_pepxml_xml(raw_dir)

        if dry_run:
            continue

        # 2) TE-only pepXML split (inside container, python+lxml)
        pepxml_files = sorted([p for p in raw_dir.iterdir() if p.is_file() and p.name.endswith(".pep.xml")])
        if not pepxml_files:
            raise FileNotFoundError(f"No pep.xml found under {raw_dir}")

        _split_pepxml_te_only_in_container(
            runner=engine_runner,
            commands_sh=outputs.commands_sh,
            log_file=outputs.log_file,
            te_fasta=te_fasta,
        )

        te_pepxml_files = sorted([p for p in raw_dir.iterdir() if p.is_file() and p.name.endswith("_TE.pep.xml")])
        if not te_pepxml_files:
            raise FileNotFoundError(f"No *_TE.pep.xml produced under {raw_dir}")

        # 3) xinteract merge + PeptideProphet (TE only)
        merged_te_pepxml = te_dir / f"{engine}_merged_TE.pep.xml"
        x_cmd = (
            f"/opt/tpp/bin/xinteract -OAPNd -PPM -eN -p0 -THREADS=16 -drev_ "
            f"-D{combined_target_decoy} "
            f"-N{merged_te_pepxml} "
            f"{raw_dir}/*_TE.pep.xml"
        )
        inner_bash = ["bash", "-lc", x_cmd]
        cmd_line = engine_runner.cmd_to_shell(engine_runner.build_run_cmd(inner_bash))
        _write_cmd(outputs.commands_sh, cmd_line)
        rc = engine_runner.run(inner_bash, dry_run=False, log_file=outputs.log_file)
        if rc != 0:
            # 139 usually means PepXMLViewer.cgi segfault; core pepXML may still be generated
            if rc == 139 and _nonempty(merged_te_pepxml):
                print(f"[WARN] xinteract exited 139 but merged pepXML exists: {merged_te_pepxml}; continuing.")
            else:
                raise RuntimeError(f"xinteract failed for {engine} (exit {rc})")

        # 4) pepxml2csv
        inner2 = ["python", "/opt/tpp/bin/pepxml2csv.py", str(merged_te_pepxml)]
        cmd_line2 = engine_runner.cmd_to_shell(engine_runner.build_run_cmd(inner2))
        _write_cmd(outputs.commands_sh, cmd_line2)
        rc2 = engine_runner.run(inner2, dry_run=False, log_file=outputs.log_file)
        if rc2 != 0:
            raise RuntimeError(f"pepxml2csv failed for {engine} (exit {rc2})")

        merged_te_csv = Path(str(merged_te_pepxml).replace(".pep.xml", ".pep.csv"))
        if not merged_te_csv.exists():
            raise FileNotFoundError(f"pepxml2csv output not found for {engine}: {merged_te_csv}")

        # 5) TE-only PSM FDR
        te_fdr_csv = _filter_by_fdr(merged_te_csv, fdr_threshold=fdr)

        if (not te_fdr_csv.exists()) or te_fdr_csv.stat().st_size == 0:
            # No TE PSM passed FDR; write an empty final marker and skip blastp
            final_marker.write_text("", encoding="utf-8")
            continue

        # 6) TE-only blastp removal vs canonical proteome
        te_filtered = _blastp_remove_canonical(
            runner=engine_runner,
            commands_sh=outputs.commands_sh,
            log_file=outputs.log_file,
            te_fdr_csv=te_fdr_csv,
            human_blastdb_prefix=human_blastdb_prefix,
            blastp_bin=blastp_bin,
            task=blastp_task,
            evalue=blastp_evalue,
            outfmt=blastp_outfmt,
            threads=blastp_threads,
        )

        shutil.copy2(te_filtered, final_marker)
