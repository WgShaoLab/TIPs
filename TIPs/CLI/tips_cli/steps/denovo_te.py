from __future__ import annotations
import ast
import math
import re
import shlex
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple
import numpy as np
import pandas as pd
import yaml

from tips_cli.runner import DockerRunner


@dataclass
class DenovoTeOutputs:
    sample_path: Path
    denovo_dir: Path
    tags_dir: Path
    te_dir: Path
    log_file: Path
    commands_sh: Path


def _nonempty(p: Path) -> bool:
    return p.exists() and p.is_file() and p.stat().st_size > 0


def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def extract_subsequences_sub(sequence: str, scores: List[float], min_score: float = 0.75, min_length: int = 8) -> List[str]:
    result: List[str] = []
    n = len(sequence)
    i = 0

    while i < n:
        if scores[i] > min_score:
            start = i
            while i < n and scores[i] > min_score:
                i += 1
            end = i

            if end - start >= min_length:
                subseq = sequence[start:end]
                subseq_scores = scores[start:end]
                result.append(subseq)

                best_avg_score = -1.0
                best_subseq = ""

                for j in range(len(subseq) - min_length + 1):
                    for k in range(j + min_length, len(subseq) + 1):
                        temp_subseq = subseq[j:k]
                        temp_scores = subseq_scores[j:k]
                        avg_score = sum(temp_scores) / len(temp_scores)
                        if avg_score > best_avg_score:
                            best_avg_score = avg_score
                            best_subseq = temp_subseq

                if best_subseq:
                    result.append(best_subseq)
        else:
            i += 1

    return result


def _sanitize_peptide(seq: str, strip_non_aa: bool = True) -> str:
    s = (seq or "").strip()
    if strip_non_aa:
        s = re.sub(r"[^A-Z]", "", s.upper())
    return s


def _parse_pepnet_file(path: Path, min_score: float, min_length: int, strip_non_aa: bool) -> List[str]:
    """
    PepNet: tab-delimited with columns including 'DENOVO' and 'Positional Score'.
    """
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return []

    if "DENOVO" not in df.columns or "Positional Score" not in df.columns:
        # Try fallback for slight header variations
        raise ValueError(f"PepNet file missing required columns: {path}")

    out: List[str] = []
    for _, row in df.iterrows():
        seq = _sanitize_peptide(str(row["DENOVO"]), strip_non_aa=strip_non_aa)
        if not seq or len(seq) < min_length:
            continue

        raw = row["Positional Score"]
        try:
            scores = ast.literal_eval(raw) if isinstance(raw, str) else raw
        except Exception:
            # Fallback to original cleaning behavior
            cleaned = re.sub(r"[\[\]\s]", "", str(raw))
            scores = [float(x) for x in cleaned.split(",") if x != ""]

        if not isinstance(scores, list) or not scores:
            continue

        subseqs = extract_subsequences_sub(seq, [float(x) for x in scores], min_score=min_score, min_length=min_length)
        out.extend(subseqs)

    return out


def _parse_instanovo_file(path: Path, min_score: float, min_length: int, strip_non_aa: bool) -> List[str]:
    """
    InstaNovo CSV:
      - peptide sequence column: preds
      - positional scores: token_log_probs (log-prob list), convert using exp() like original code.
    """
    df = pd.read_csv(path)
    if df.empty:
        return []

    if "preds" not in df.columns or "token_log_probs" not in df.columns:
        raise ValueError(f"InstaNovo file missing required columns (preds, token_log_probs): {path}")

    out: List[str] = []
    for _, row in df.iterrows():
        seq = _sanitize_peptide(str(row["preds"]), strip_non_aa=strip_non_aa)
        if not seq or len(seq) < min_length:
            continue

        raw = row["token_log_probs"]
        try:
            logps = ast.literal_eval(raw) if isinstance(raw, str) else raw
        except Exception:
            # Fallback to original slicing behavior: x[1:-1].split(',')
            s = str(raw).strip()
            if s.startswith("[") and s.endswith("]"):
                s = s[1:-1]
            logps = [float(x) for x in s.split(",") if x.strip() != ""]

        if not isinstance(logps, list) or not logps:
            continue

        scores = [float(math.exp(float(v))) for v in logps]
        subseqs = extract_subsequences_sub(seq, scores, min_score=min_score, min_length=min_length)
        out.extend(subseqs)

    return out


def _read_casanovo_mztab_without_mtd(path: Path) -> pd.DataFrame:
    """
    Match the original behavior: remove lines starting with 'MTD', then pandas read_csv with tab separator.
    """
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines(True)
    filtered_lines = [ln for ln in lines if not ln.startswith("MTD")]
    if not filtered_lines:
        return pd.DataFrame()
    from io import StringIO
    return pd.read_csv(StringIO("".join(filtered_lines)), sep="\t")


def _parse_casanovo_mztab(path: Path, min_score: float, min_length: int, strip_non_aa: bool) -> List[str]:
    """
    Casanovo mztab:
      - sequence column: sequence
      - AA scores column: opt_ms_run[1]_aa_scores (comma-separated floats)
      - use extract_subsequences_sub identical to original.
    """
    df = _read_casanovo_mztab_without_mtd(path)
    if df.empty:
        return []

    if "sequence" not in df.columns or "opt_ms_run[1]_aa_scores" not in df.columns:
        raise ValueError(f"Casanovo mztab missing required columns: {path}")

    out: List[str] = []
    for _, row in df.iterrows():
        seq = _sanitize_peptide(str(row["sequence"]), strip_non_aa=strip_non_aa)
        if not seq or len(seq) < min_length:
            continue

        raw = row["opt_ms_run[1]_aa_scores"]
        if pd.isna(raw):
            continue
        try:
            scores = [float(x) for x in str(raw).split(",") if str(x).strip() != ""]
        except Exception:
            continue

        subseqs = extract_subsequences_sub(seq, scores, min_score=min_score, min_length=min_length)
        out.extend(subseqs)

    return out


def _write_fasta(peptides: Iterable[str], out_fa: Path) -> None:
    _ensure_dir(out_fa.parent)
    with out_fa.open("w", encoding="utf-8") as f:
        for pep in peptides:
            p = str(pep).strip()
            if not p:
                continue
            f.write(f">{p}\n{p}\n")


def _read_fasta_headers_as_peptides(fa: Path) -> List[str]:
    txt = fa.read_text(encoding="utf-8", errors="ignore")
    parts = txt.split(">")[1:]
    out: List[str] = []
    for rec in parts:
        head = rec.splitlines()[0].strip()
        if head:
            out.append(head)
    return out


def _write_split_fastas(query_fa: Path, out_dir: Path, num_files: int) -> List[Path]:
    peptides = _read_fasta_headers_as_peptides(query_fa)
    if not peptides:
        return []

    num_files = max(1, int(num_files))
    chunk = len(peptides) // num_files
    rem = len(peptides) % num_files

    out_paths: List[Path] = []
    start = 0
    for i in range(num_files):
        end = start + chunk + (1 if i < rem else 0)
        sub = peptides[start:end]
        start = end
        if not sub:
            continue
        p = out_dir / f"sub_fasta_{i + 1}.fasta"
        _write_fasta(sub, p)
        out_paths.append(p)

    return out_paths


def _append_cmd(commands_sh: Path, cmd_list: List[str]) -> None:
    commands_sh.parent.mkdir(parents=True, exist_ok=True)
    line = " ".join(shlex.quote(x) for x in cmd_list)
    with commands_sh.open("a", encoding="utf-8") as f:
        f.write(line + "\n")


def _run_one_blastp(
    runner: DockerRunner,
    blastp_bin: str,
    db_prefix: str,
    query_fa: Path,
    out_txt: Path,
    outfmt: str,
    evalue: float,
    threads: int,
    dry_run: bool,
    commands_sh: Path,
    log_file: Optional[Path],
) -> int:
    """
    Run blastp-short inside docker (tips-core).
    """
    inner_cmd = [
        blastp_bin,
        "-task",
        "blastp-short",
        "-query",
        str(query_fa),
        "-db",
        str(db_prefix),
        "-out",
        str(out_txt),
        "-outfmt",
        outfmt,
        "-evalue",
        str(evalue),
        "-num_threads",
        str(int(threads)),
    ]

    full = runner.build_run_cmd(inner_cmd)
    _append_cmd(commands_sh, full)

    return runner.run(inner_cmd, dry_run=dry_run, log_file=log_file)


from pathlib import Path
import pandas as pd


def _blastp_merge_and_filter(child_dir: Path, min_pident: float) -> Path:
    """
    Merge sub blast outputs and generate:
      - merged_blastp.txt
      - merged_blastp_filter.txt (coarse filter: qstart==1, gapopen==0, pident>=min_pident)
      - merged_blastp_TE.csv (legacy logic with L_count, mismatch<=L_count, and I/L identity by I_L_identity)
    """
    txts = sorted([p for p in child_dir.glob("sub_fasta_*.txt") if p.is_file()])
    merged_txt = child_dir / "merged_blastp.txt"
    merged_flt = child_dir / "merged_blastp_filter.txt"
    merged_csv = child_dir / "merged_blastp_TE.csv"

    # 1) Merge raw blast outputs (cat-like, guarantee newline between files)
    with merged_txt.open("wb") as out:
        for fp in txts:
            with fp.open("rb") as fin:
                out.write(fin.read())
            out.write(b"\n")

    if not merged_txt.exists() or merged_txt.stat().st_size == 0:
        pd.DataFrame().to_csv(merged_csv, index=False, sep="\t")
        return merged_csv

    # 2) Coarse filter (awk-like: split by ANY whitespace)
    min_pident_f = float(min_pident)

    with merged_txt.open("r", encoding="utf-8", errors="ignore") as fin, merged_flt.open(
        "w", encoding="utf-8"
    ) as fout:
        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue

            # Mimic awk default field splitting (any whitespace)
            parts = line.split()
            if len(parts) < 14:
                continue

            try:
                pident = float(parts[2])
                gapopen = int(float(parts[5]))
                qstart = int(float(parts[6]))
            except ValueError:
                continue

            if qstart == 1 and gapopen == 0 and pident >= min_pident_f:
                # Write a consistent, tab-delimited line for downstream pandas
                fout.write("\t".join(parts[:14]) + "\n")

    if not merged_flt.exists() or merged_flt.stat().st_size == 0:
        pd.DataFrame().to_csv(merged_csv, index=False, sep="\t")
        return merged_csv

    # 3) Generate merged_blastp_TE.csv (legacy pandas logic)
    def I_L_identity(q_seq: str, t_seq: str) -> bool:
        """Legacy I/L identity logic from the historical implementation."""
        if q_seq == t_seq:
            return True

        diff_list = []
        for i in range(len(q_seq)):
            if q_seq[i] != t_seq[i] and set([q_seq[i], t_seq[i]]) == set(["I", "L"]):
                diff_list.append(t_seq[i])

        diff_list = list(set(diff_list))
        if not diff_list:
            return False

        return "".join(sorted(diff_list)) in "IL"

    column_names = [
        "qaccver",
        "saccver",
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
        "qseq",
        "sseq",
    ]

    # Prefer tab-delimited; if it collapses to a single column, fall back to whitespace parsing
    df = pd.read_csv(merged_flt, names=column_names, sep="\t", header=None)
    if df.shape[1] == 1:
        df = pd.read_csv(merged_flt, names=column_names, sep=r"\s+", engine="python", header=None)

    if df.empty:
        df.to_csv(merged_csv, index=False, sep="\t")
        return merged_csv

    # Ensure numeric columns are numeric to make comparisons robust
    for c in ["pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # Legacy helper columns
    df["qaccver"] = df["qaccver"].astype(str)
    df["qseq"] = df["qseq"].astype(str)
    df["sseq"] = df["sseq"].astype(str)

    df["q_length"] = df["qaccver"].str.len()
    df["L_count"] = df["qaccver"].str.count("L")

    # Legacy filters
    df_filter = df[df["length"] == df["q_length"]]
    df_filter = df_filter[df_filter["qseq"].str.len() == df_filter["sseq"].str.len()]
    df_filter = df_filter[df_filter["gapopen"] == 0]
    df_filter = df_filter[df_filter["mismatch"] <= df_filter["L_count"]]

    # Legacy I/L identity definition
    df_filter["I/L_identity"] = df_filter.apply(
        lambda x: I_L_identity(x["qseq"], x["sseq"]), axis=1
    )

    # Legacy output schema (no qseq_I2L/sseq_I2L; include L_count)
    out_cols = [
        "qaccver",
        "saccver",
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
        "qseq",
        "sseq",
        "q_length",
        "L_count",
        "I/L_identity",
    ]
    df_filter.to_csv(merged_csv, index=False, sep="\t", columns=out_cols)
    return merged_csv

def denovo_te_step(
    outputs: DenovoTeOutputs,
    image: str,
    config_path: Path,
    tools: Sequence[str],
    blastp_workers: int,
    threads_per_task: int,
    dry_run: bool,
    resume: bool,
    overwrite: bool,
    keep_intermediate: bool,
) -> None:
    sample_path = outputs.sample_path.resolve()
    denovo_dir = outputs.denovo_dir.resolve()
    _ensure_dir(outputs.tags_dir)
    _ensure_dir(outputs.te_dir)
    outputs.commands_sh.write_text("#!/usr/bin/env bash\nset -euo pipefail\n\n", encoding="utf-8")

    if not config_path.is_file():
        raise FileNotFoundError(f"Config not found: {config_path}")

    cfg = yaml.safe_load(config_path.read_text(encoding="utf-8"))

    import os
    from pathlib import Path

    def _abspath_no_realpath(p: Path) -> Path:
        # os.path.abspath does NOT dereference symlinks (unlike Path.resolve)
        return Path(os.path.abspath(str(p)))

    def _resolve_cfg_path(v: str, sample_path: Path, config_path: Path) -> Path:
        p = Path(v).expanduser()
        if p.is_absolute():
            return p

        # Prefer sample-relative path to make each sample self-contained and reproducible.
        p_sample = _abspath_no_realpath(sample_path / p)
        if p_sample.exists():
            return p_sample

        # Fallback to config-relative path for shared configs.
        p_cfg = _abspath_no_realpath(config_path.parent / p)
        if p_cfg.exists():
            return p_cfg

        # If neither exists, return sample-relative so downstream errors point to the sample workspace.
        return p_sample

    te_protein_npy = _resolve_cfg_path(str(cfg["te"]["protein_npy"]), sample_path, config_path)
    te_class_npy = _resolve_cfg_path(str(cfg["te"]["class_npy"]), sample_path, config_path)
    db_prefix_path = _resolve_cfg_path(str(cfg["te"]["blast_db_prefix"]), sample_path, config_path)
    db_prefix = str(db_prefix_path)

    def _blastdb_exists(prefix: Path) -> bool:
        parent = prefix.parent
        base = prefix.name
        if not parent.exists():
            return False

        # Single-volume pattern
        if (parent / f"{base}.pin").exists() or (parent / f"{base}.phr").exists() or (parent / f"{base}.psq").exists():
            return True

        # Multi-volume pattern: <prefix>.<NN>.pin/.phr/.psq
        # Accept any numeric NN (00,01,02,...)
        return bool(
            list(parent.glob(f"{base}.[0-9][0-9].pin"))
            or list(parent.glob(f"{base}.[0-9][0-9].phr"))
            or list(parent.glob(f"{base}.[0-9][0-9].psq"))
        )

    if not _blastdb_exists(db_prefix_path):
        raise FileNotFoundError(
            f"BLAST database not found for prefix: {db_prefix_path}. "
            f"Expected files like {db_prefix_path}.pin/.phr/.psq or {db_prefix_path}.00.pin/.00.phr/.00.psq"
        )
    db_prefix = str(db_prefix_path)

    blastp_cfg = cfg.get("blastp", {})
    blastp_bin = str(blastp_cfg.get("bin", "blastp"))
    outfmt = str(
        blastp_cfg.get(
            "outfmt",
            "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq",
        )
    )
    evalue = float(blastp_cfg.get("evalue", 20000))
    min_pident = float(blastp_cfg.get("min_pident", 75))

    tag_cfg = cfg.get("tag_extract", {})
    min_score = float(tag_cfg.get("min_score", 0.75))
    min_length = int(tag_cfg.get("min_length", 8))
    strip_non_aa = bool(tag_cfg.get("strip_non_aa", True))

    tool_globs = cfg.get("tools", {})

    # Outputs (fixed naming with neutral label: Merged)
    query_fa = outputs.tags_dir / "Denovo_result_merged_Merged.fasta"
    child_dir = outputs.tags_dir / "Denovo_result_merged_Merged_Child_fasta"
    merged_blast_csv = child_dir / "merged_blastp_TE.csv"
    out_te_fa = outputs.te_dir / "Denovo_TE_SoftMerged_Merged.fasta"
    out_te_tsv = outputs.te_dir / "Denovo_TE_SoftMerged_Merged.tsv"

    if overwrite:
        for p in [query_fa, out_te_fa, out_te_tsv]:
            if p.exists():
                p.unlink()

    # 1) Extract per-tool tags from raw result files and write merged_raw_peptide.fasta (for transparency)
    peptide_soft: Dict[str, List[str]] = {}

    for tool in tools:
        tool_l = tool.strip().lower()
        if not tool_l:
            continue

        in_dir = denovo_dir / tool_l
        if not in_dir.exists():
            continue

        glob_pat = "*"
        if isinstance(tool_globs, dict) and tool_l in tool_globs and isinstance(tool_globs[tool_l], dict):
            glob_pat = str(tool_globs[tool_l].get("glob", "*"))

        files = sorted([p for p in in_dir.glob(glob_pat) if p.is_file()])
        if not files:
            continue

        peptides: List[str] = []
        for fp in files:
            if tool_l == "pepnet":
                peptides.extend(_parse_pepnet_file(fp, min_score=min_score, min_length=min_length, strip_non_aa=strip_non_aa))
            elif tool_l == "instanovo":
                peptides.extend(
                    _parse_instanovo_file(fp, min_score=min_score, min_length=min_length, strip_non_aa=strip_non_aa)
                )
            elif tool_l == "casanovo":
                peptides.extend(
                    _parse_casanovo_mztab(fp, min_score=min_score, min_length=min_length, strip_non_aa=strip_non_aa)
                )
            else:
                raise ValueError(f"Unsupported tool for denovo-te: {tool}")

        peptides = [p for p in peptides if p and len(p) >= min_length]
        peptides = list(set(peptides))

        per_tool_fa = outputs.tags_dir / tool_l / "merged_raw_peptide.fasta"
        if overwrite and per_tool_fa.exists():
            per_tool_fa.unlink()
        if not (resume and _nonempty(per_tool_fa)):
            _write_fasta(peptides, per_tool_fa)

        for pep in peptides:
            peptide_soft.setdefault(pep, []).append(tool_l)

    if not peptide_soft:
        raise FileNotFoundError(f"No de novo peptides extracted under: {denovo_dir}")

    # 2) Write merged query fasta (header=peptide, seq=peptide)
    if not (resume and _nonempty(query_fa)):
        _write_fasta(peptide_soft.keys(), query_fa)

    # 3) Run blastp-short in docker (split + parallel)
    import os

    def _real_parent_dir(p: Path) -> Path:
        # realpath follows symlinks -> points to the true location on host
        rp = Path(os.path.realpath(str(p)))
        return rp.parent

    extra_mounts = set()

    # mount real directories of TE resources (npy + blastdb prefix directory)
    extra_mounts.add(_real_parent_dir(te_protein_npy))
    extra_mounts.add(_real_parent_dir(te_class_npy))
    extra_mounts.add(_real_parent_dir(db_prefix_path))

    # Also mount the real directory of any blastdb volume files if prefix itself is a symlink file
    # (safe even if not)
    # Note: db_prefix_path may not exist as a file; it's a prefix. Mount its directory anyway:
    extra_mounts.add(Path(os.path.realpath(str(db_prefix_path.parent))))

    runner = DockerRunner(
        image=image,
        mount_root=sample_path,
        user_map=True,
        workdir_in_container="/tmp",
        extra_mounts=sorted(extra_mounts),
    )

    if not (resume and _nonempty(merged_blast_csv)):
        if child_dir.exists() and overwrite:
            subprocess.run(["rm", "-rf", str(child_dir)], check=False)
        child_dir.mkdir(parents=True, exist_ok=True)

        sub_fastas = _write_split_fastas(query_fa, child_dir, num_files=max(1, blastp_workers))

        futures = []
        with ThreadPoolExecutor(max_workers=max(1, int(blastp_workers))) as ex:
            for sub_fa in sub_fastas:
                out_txt = sub_fa.with_suffix(".txt")
                futures.append(
                    ex.submit(
                        _run_one_blastp,
                        runner,
                        blastp_bin,
                        db_prefix,
                        sub_fa,
                        out_txt,
                        outfmt,
                        evalue,
                        threads_per_task,
                        dry_run,
                        outputs.commands_sh,
                        outputs.log_file,
                    )
                )
            for fut in as_completed(futures):
                rc = fut.result()
                if rc != 0:
                    raise RuntimeError(f"blastp failed (rc={rc}). See: {outputs.log_file}")

        if dry_run:
            return

        _blastp_merge_and_filter(child_dir, min_pident=min_pident)

    # 4) Build TE FASTA
    if not te_protein_npy.is_file():
        raise FileNotFoundError(f"TE protein npy not found: {te_protein_npy}")
    if not te_class_npy.is_file():
        raise FileNotFoundError(f"TE class npy not found: {te_class_npy}")

    te_pro_dic = np.load(str(te_protein_npy), allow_pickle=True).item()
    te_class_dic = np.load(str(te_class_npy), allow_pickle=True).item()

    df = pd.read_csv(merged_blast_csv, sep="\t")
    if df.empty:
        out_te_fa.write_text("", encoding="utf-8")
        out_te_tsv.write_text("", encoding="utf-8")
        return

    if "I/L_identity" in df.columns:
        df = df[df["I/L_identity"] == True]  # noqa: E712

    df["soft"] = df["qaccver"].astype(str).apply(lambda x: peptide_soft.get(x, []))

    def _te_name_from_saccver(x: str) -> str:
        s = str(x)
        if "Class=DNA;" not in s:
            return s.split("::")[0]
        return s.split(";")[0].split("=")[1]

    df["TE_name"] = df["saccver"].astype(str).apply(_te_name_from_saccver)
    df["TE_class"] = df["TE_name"].astype(str).apply(lambda x: te_class_dic[x] if x in te_class_dic else "DNA")

    df = df.drop_duplicates(subset=["qaccver", "saccver"])

    out_df = df[["saccver", "TE_class", "soft", "TE_name"]].copy()
    out_df = (
        out_df.groupby(["saccver", "TE_class", "TE_name"], as_index=False)
        .agg({"soft": lambda xs: sum(xs, [])})
    )
    out_df["soft"] = out_df["soft"].apply(lambda xs: "_".join(sorted(set(xs))))

    _ensure_dir(out_te_fa.parent)
    out_df.to_csv(out_te_tsv, index=False, sep="\t")

    with out_te_fa.open("w", encoding="utf-8") as f:
        for _, row in out_df.iterrows():
            soft = row["soft"]
            te_header = row["saccver"]
            te_name = row["TE_name"]
            te_class = row["TE_class"]
            if te_header not in te_pro_dic:
                continue
            seq = te_pro_dic[te_header]
            header = f"Denovo|{soft}:::{te_name}:::{te_class}:::{te_header}"
            f.write(f">{header}\n{seq}\n")

    if not keep_intermediate:
        subprocess.run(["rm", "-rf", str(child_dir)], check=False)
