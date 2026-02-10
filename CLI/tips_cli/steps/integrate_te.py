from __future__ import annotations

import shlex
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Set, List

import pandas as pd

from tips_cli.runner import DockerRunner


@dataclass
class IntegrateTeOutputs:
    sample_path: Path
    integrate_dir: Path
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


def _write_text(p: Path, text: str) -> None:
    _ensure_dir(p.parent)
    p.write_text(text, encoding="utf-8")


def _peptide_fasta_from_tsv(
    peptide_tsv: Path,
    out_fasta: Path,
    peptide_col_candidates: List[str],
    length_lt: Optional[int] = None,
) -> List[str]:
    """
    Convert a Philosopher peptide.tsv into a FASTA file.

    Returns a sorted unique peptide list used to write the FASTA.
    """
    df = pd.read_csv(peptide_tsv, sep="\t")
    pep_col = None
    for c in peptide_col_candidates:
        if c in df.columns:
            pep_col = c
            break
    if pep_col is None:
        raise ValueError(
            f"Cannot find peptide column in {peptide_tsv}. Tried: {peptide_col_candidates}"
        )

    peptides: Set[str] = set()
    for x in df[pep_col].dropna().astype(str).tolist():
        s = x.strip()
        if not s:
            continue
        if length_lt is not None and not (0 < len(s) < length_lt):
            continue
        peptides.add(s)

    pep_list = sorted(peptides)

    lines: List[str] = []
    for p in pep_list:
        lines.append(f">{p}\n{p}\n")
    _write_text(out_fasta, "".join(lines))
    return pep_list


def _strict_blastp_remove_set(blast_outfmt6: Path) -> Set[str]:
    """
    Parse BLAST outfmt 6 and return peptides to remove using strict rules:
      - pident == 100
      - qstart == 1
      - qend == query length
    """
    if not _nonempty(blast_outfmt6):
        return set()

    cols = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
    ]
    bdf = pd.read_csv(blast_outfmt6, sep="\t", names=cols)
    bdf["qseqid"] = bdf["qseqid"].astype(str)
    bdf["qlen"] = bdf["qseqid"].apply(len)

    strict = bdf[
        (bdf["pident"] == 100.0)
        & (bdf["qstart"] == 1)
        & (bdf["qend"] == bdf["qlen"])
    ]
    return set(strict["qseqid"].drop_duplicates().tolist())


def integrate_te_step(
    outputs: IntegrateTeOutputs,
    image: str,
    combined_target_decoy_fasta: Path,
    pepxml_comet: Path,
    pepxml_msfragger: Path,
    pepxml_msgfplus: Path,
    fdr: float = 0.03,
    threads: int = 20,
    blastp_db_prefix: str = "resources/db/human_blastpdb/human_proteome_UP000005640_noERV_LINE",
    blastp_threads: int = 20,
    hla_alleles: Optional[str] = None,
    dry_run: bool = False,
    resume: bool = False,
    overwrite: bool = False,
) -> None:
    """
    Integrate TE-only pepXML from multiple search engines.

    Pipeline (aligned with SoftPipeline_class_iProphet.py):
      1) philosopher workspace --init
      2) philosopher database --annotate <combined_target_decoy_fasta>
      3) RefreshParser for Comet/MSFragger pepXML
      4) philosopher iprophet (3 pepXML) -> proteinprophet -> filter (--pep/--psm) -> report
      5) blastp strict full-length identity to remove canonical peptides
      6) optional MixMHCpred binder prediction (if HLA alleles provided)

    All outputs are written under:
      <sample_path>/integrate/
    """
    sample_path = outputs.sample_path.resolve()
    integrate_dir = outputs.integrate_dir.resolve()

    # Flat layout: keep everything directly under integrate_dir (TE-only mode).
    work_dir = integrate_dir
    _ensure_dir(work_dir)

    in_dir = work_dir
    ph_dir = work_dir
    bp_dir = work_dir
    bind_dir = work_dir

    final_table = work_dir / "final_TE_results.tsv"

    if resume and _nonempty(final_table):
        return

    if overwrite and integrate_dir.exists():
        # Remove only known integrate outputs; never touch search outputs.
        patterns = [
            "comet_merged_TE.pep.xml",
            "msfragger_merged_TE.pep.xml",
            "msgfplus_merged_TE.pep.xml",
            "interact.*",
            "*.pep.xml.index",
            "peptide.tsv",
            "psm.tsv",
            "protein.tsv",
            "ion.tsv",
            "peptide_blastp.fasta",
            "peptide_blastp.outfmt6.tsv",
            "binding.fasta",
            "mixmhcpred.out.tsv",
            "final_TE_results.tsv",
            ".meta",
        ]
        for pat in patterns:
            for p in integrate_dir.glob(pat):
                try:
                    if p.is_dir():
                        continue
                    p.unlink()
                except Exception:
                    pass


    # Link inputs into integrate folder for stable, reproducible paths
    pep1 = in_dir / "comet_merged_TE.pep.xml"
    pep2 = in_dir / "msfragger_merged_TE.pep.xml"
    pep3 = in_dir / "msgfplus_merged_TE.pep.xml"
    _safe_symlink(pepxml_comet.resolve(), pep1)
    _safe_symlink(pepxml_msfragger.resolve(), pep2)
    _safe_symlink(pepxml_msgfplus.resolve(), pep3)

    # Resolve BLAST DB prefix to an absolute path inside the mounted sample directory.
    # BLAST expects index/alias files next to the prefix (e.g., .pin/.psq/.phr or .pdb/.pot/.ptf/.pto).
    db_prefix_path = Path(blastp_db_prefix)
    if not db_prefix_path.is_absolute():
        db_prefix_path = (sample_path / db_prefix_path).resolve()
    else:
        db_prefix_path = db_prefix_path.resolve()


    # If BLAST DB index files are symlinks pointing outside sample_path, mount the real DB directory.
    pin = db_prefix_path.parent / f"{db_prefix_path.name}.pin"
    psq = db_prefix_path.parent / f"{db_prefix_path.name}.psq"
    phr = db_prefix_path.parent / f"{db_prefix_path.name}.phr"
    idx_files = [p for p in (pin, psq, phr) if p.exists()]
    if len(idx_files) == 0:
        raise FileNotFoundError(
            f"BLAST DB index files not found for prefix: {db_prefix_path} "
            f"(expected .pin/.psq/.phr next to the prefix)."
        )

    # Resolve the real directory where the BLAST DB files actually live
    real_db_dir = idx_files[0].resolve().parent
    real_db_prefix = real_db_dir / db_prefix_path.name

    # DockerRunner only supports same-path mounts (-v host:host), so mount the real DB dir as-is
    extra_mounts: List[Path] = []
    try:
        real_db_dir.relative_to(sample_path)
        # Real DB already under sample_path -> visible via mount_root
    except ValueError:
        extra_mounts.append(real_db_dir)

    # Use the real DB prefix path inside container (same as host path due to same-path mount)
    db_prefix_in_container = str(real_db_prefix)


    runner = DockerRunner(
        image=image,
        mount_root=sample_path,
        user_map=True,
        workdir_in_container=str(sample_path),
        extra_mounts=extra_mounts,
    )


    def _run(inner_bash: List[str], step_name: str) -> None:
        cmd_line = runner.cmd_to_shell(runner.build_run_cmd(inner_bash))
        _write_cmd(outputs.commands_sh, cmd_line)
        rc = runner.run(inner_bash, dry_run=dry_run, log_file=outputs.log_file)
        if rc != 0:
            raise RuntimeError(f"integrate-te {step_name} failed (exit {rc})")

    # 1) Philosopher workspace init + DB annotate
    _run(
        [
            "bash",
            "-lc",
            f"""
set -euo pipefail
cd {shlex.quote(str(ph_dir))}
philosopher workspace --init
philosopher database --annotate {shlex.quote(str(combined_target_decoy_fasta))}
""",
        ],
        "philosopher-init",
    )

    # 2) RefreshParser (Comet + MSFragger only)
    _run(
        [
            "bash",
            "-lc",
            f"""
set -euo pipefail
cd {shlex.quote(str(ph_dir))}
RefreshParser {shlex.quote(str(pep1))} {shlex.quote(str(combined_target_decoy_fasta))}
RefreshParser {shlex.quote(str(pep2))} {shlex.quote(str(combined_target_decoy_fasta))}
""",
        ],
        "refreshparser",
    )

    # 3) iProphet + proteinprophet + filter + report
    _run(
        [
            "bash",
            "-lc",
            f"""
set -euo pipefail
cd {shlex.quote(str(ph_dir))}
philosopher iprophet --threads {int(threads)} {shlex.quote(str(pep1))} {shlex.quote(str(pep2))} {shlex.quote(str(pep3))}
philosopher proteinprophet --iprophet interact.iproph.pep.xml
philosopher filter --pepxml interact.iproph.pep.xml --pep {float(fdr)} --psm {float(fdr)}
philosopher report
""",
        ],
        "philosopher-iprophet",
    )

    peptide_tsv = ph_dir / "peptide.tsv"
    if not peptide_tsv.exists():
        raise FileNotFoundError(f"Missing Philosopher output: {peptide_tsv}")

    def _filter_sp_mapped_peptides(
            peptide_tsv: Path,
            out_tsv: Path,
            mapped_protein_col_candidates: List[str] = None,
    ) -> None:
        """
        Remove peptides whose mapped proteins contain 'sp|' (canonical UniProt entries).

        This filter is applied BEFORE blastp and binder prediction to eliminate
        obvious canonical-mapped peptides early.
        """
        if mapped_protein_col_candidates is None:
            mapped_protein_col_candidates = [
                "Mapped Proteins",
                "Mapped_Proteins",
                "Proteins",
                "Protein",
            ]

        df = pd.read_csv(peptide_tsv, sep="\t")

        mp_col = None
        for c in mapped_protein_col_candidates:
            if c in df.columns:
                mp_col = c
                break

        # If no mapped-protein column exists, keep all peptides (fail-safe)
        if mp_col is None:
            df.to_csv(out_tsv, sep="\t", index=False)
            return

        # Remove rows where mapped proteins contain canonical UniProt entries
        mask_sp = df[mp_col].astype(str).str.contains(r"sp\|", regex=True)
        df_f = df[~mask_sp].copy()

        df_f.to_csv(out_tsv, sep="\t", index=False)

    # Pre-filter peptides mapped to canonical UniProt proteins (sp|)
    peptide_prefilter_tsv = ph_dir / "peptide.prefilter.tsv"
    _filter_sp_mapped_peptides(
        peptide_tsv=peptide_tsv,
        out_tsv=peptide_prefilter_tsv,
    )

    # Use prefiltered peptide table for downstream blastp and binder
    peptide_tsv = peptide_prefilter_tsv

    # 4) blastp strict canonical removal
    blast_fa = bp_dir / "peptide_blastp.fasta"
    blast_out = bp_dir / "peptide_blastp.outfmt6.tsv"
    _peptide_fasta_from_tsv(
        peptide_tsv=peptide_tsv,
        out_fasta=blast_fa,
        peptide_col_candidates=["peptide", "Peptide"],
        length_lt=None,
    )

    _run(
        [
            "bash",
            "-lc",
            f"""
set -euo pipefail
cd {shlex.quote(str(bp_dir))}
blastp -task blastp-short \
  -query {shlex.quote(str(blast_fa))} \
  -db {shlex.quote(str(db_prefix_in_container))} \
  -out {shlex.quote(str(blast_out))} \
  -outfmt 6 -evalue 20000 -num_threads {int(blastp_threads)}
""",
        ],
        "blastp",
    )

    remove_set = _strict_blastp_remove_set(blast_out)

    df = pd.read_csv(peptide_tsv, sep="\t")
    pep_col = "peptide" if "peptide" in df.columns else ("Peptide" if "Peptide" in df.columns else None)
    if pep_col is None:
        raise ValueError("peptide.tsv missing peptide column")

    df_f = df[~df[pep_col].astype(str).isin(remove_set)].copy()

    # 5) optional MixMHCpred binder prediction
    if hla_alleles:
        bind_fa = bind_dir / "binding.fasta"
        bind_out = bind_dir / "mixmhcpred.out.tsv"

        # MixMHCpred is typically used for HLA-I (8-14); follow legacy behavior: len < 15
        peptides = sorted(set([p for p in df_f[pep_col].astype(str).tolist() if 0 < len(p) < 15]))
        lines = []
        for p in peptides:
            lines.append(f">{p}\n{p}\n")
        _write_text(bind_fa, "".join(lines))

        _run(
            [
                "bash",
                "-lc",
                f"""
set -euo pipefail
cd {shlex.quote(str(bind_dir))}
MixMHCpred -i {shlex.quote(str(bind_fa))} -o {shlex.quote(str(bind_out))} -a {shlex.quote(str(hla_alleles))}
""",
            ],
            "mixmhcpred",
        )

        if _nonempty(bind_out):
            rdf = pd.read_csv(bind_out, sep="\t", comment="#")
            # Be defensive for different MixMHCpred outputs
            if {"Peptide", "%Rank_bestAllele", "BestAllele"}.issubset(set(rdf.columns)):
                sub = rdf[["Peptide", "%Rank_bestAllele", "BestAllele"]].drop_duplicates()
                out_merge = df_f.merge(sub, left_on=pep_col, right_on="Peptide", how="left")
                out_merge.drop(columns=["Peptide"], inplace=True, errors="ignore")
            else:
                out_merge = df_f
        else:
            out_merge = df_f

        out_merge.to_csv(final_table, sep="\t", index=False)
    else:
        df_f.to_csv(final_table, sep="\t", index=False)
