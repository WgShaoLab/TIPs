from __future__ import annotations

import ast
import shlex
import shutil
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


def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _nonempty(p: Path) -> bool:
    return p.exists() and p.is_file() and p.stat().st_size > 0


def _write_cmd(commands_sh: Path, cmd_line: str) -> None:
    _ensure_dir(commands_sh.parent)
    with commands_sh.open("a", encoding="utf-8") as f:
        f.write(cmd_line.rstrip() + "\n")


def _write_text(p: Path, text: str) -> None:
    _ensure_dir(p.parent)
    p.write_text(text, encoding="utf-8")


def _filter_by_fdr(pepxml_csv: Path, fdr_threshold: float) -> Path:
    """使用原生 Python 逻辑复现 TPP/Philosopher 的 PSM FDR 过滤，完全摒弃 philosopher filter"""
    out_path = pepxml_csv.with_name(pepxml_csv.name.replace(".pep.csv", "_fdr.csv"))

    try:
        df_all = pd.read_csv(pepxml_csv, sep="\t")
    except Exception:
        pd.DataFrame({"peptide": []}).to_csv(out_path, index=False, sep="\t")
        return out_path

    if df_all.empty or "protein" not in df_all.columns or "analysis_result" not in df_all.columns:
        pd.DataFrame({"peptide": []}).to_csv(out_path, index=False, sep="\t")
        return out_path

    df_all["protein_primary"] = df_all["protein"].apply(lambda x: str(x).split("#")[0])

    def _extract_prob(val):
        if pd.isna(val):
            return 0.0
        try:
            res = ast.literal_eval(str(val))
            if isinstance(res, list):
                for item in res:
                    if isinstance(item, dict):
                        # iProphet 跑完后生成的是 interprophet_result
                        if "interprophet_result" in item:
                            return float(item["interprophet_result"]["probability"])
                        # 兼容只跑了单引擎没跑 iProphet 的情况
                        if "peptideprophet_result" in item:
                            return float(item["peptideprophet_result"]["probability"])
        except Exception:
            pass
        return 0.0

    df_all["probability"] = df_all["analysis_result"].apply(_extract_prob)
    df_all["decoy"] = df_all["protein_primary"].apply(lambda x: 1 if "rev_" in str(x) else 0)
    df_all = df_all.sort_values(by=["probability"], ascending=False).reset_index(drop=True)

    df_merged = pd.DataFrame()
    if "assumed_charge" in df_all.columns:
        for charge, df_group in df_all.groupby("assumed_charge"):
            df = df_group.copy()
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
    else:
        df = df_all.copy()
        df["target_cumsum"] = (df["decoy"] == 0).cumsum()
        df["decoy_cumsum"] = (df["decoy"] == 1).cumsum()
        df = df[df["target_cumsum"] > 0].copy()
        if not df.empty:
            df["fdr"] = df["decoy_cumsum"] / df["target_cumsum"]
            df.reset_index(drop=True, inplace=True)
            reversed_df = df.iloc[::-1]
            ok = reversed_df[reversed_df["fdr"] < fdr_threshold]
            if len(ok) > 0:
                first_index = ok.index[0]
                filtered_df = df.iloc[: first_index + 1]
                df_merged = filtered_df[filtered_df["decoy"] == 0]

    if df_merged.empty:
        pd.DataFrame({"peptide": []}).to_csv(out_path, index=False, sep="\t")
        return out_path

    df_merged = df_merged.sort_values(by=["probability"], ascending=False)
    df_merged.to_csv(out_path, index=False, sep="\t")
    return out_path


def _strict_blastp_remove_set(blast_outfmt6: Path) -> Set[str]:
    if not _nonempty(blast_outfmt6):
        return set()

    cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    bdf = pd.read_csv(blast_outfmt6, sep="\t", names=cols)
    bdf["qseqid"] = bdf["qseqid"].astype(str)
    bdf["qlen"] = bdf["qseqid"].apply(len)

    strict = bdf[(bdf["pident"] == 100.0) & (bdf["qstart"] == 1) & (bdf["qend"] == bdf["qlen"])]
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
    sample_path = outputs.sample_path.resolve()
    integrate_dir = outputs.integrate_dir.resolve()

    work_dir = integrate_dir
    _ensure_dir(work_dir)

    final_table = work_dir / "final_TE_results.tsv"

    if resume and _nonempty(final_table) and not overwrite:
        return

    if overwrite and integrate_dir.exists():
        for p in integrate_dir.glob("*"):
            if p.is_file() or p.is_symlink():
                p.unlink()
        meta_dir = integrate_dir / ".meta"
        if meta_dir.exists() and meta_dir.is_dir():
            shutil.rmtree(meta_dir, ignore_errors=True)

    pep1 = work_dir / "comet_merged_TE.pep.xml"
    pep2 = work_dir / "msfragger_merged_TE.pep.xml"
    pep3 = work_dir / "msgfplus_merged_TE.pep.xml"

    if not dry_run:
        for p_in, p_out in [(pepxml_comet, pep1), (pepxml_msfragger, pep2), (pepxml_msgfplus, pep3)]:
            if p_out.exists() or p_out.is_symlink():
                p_out.unlink(missing_ok=True)
            shutil.copy2(p_in.resolve(), p_out)

    db_prefix_path = Path(blastp_db_prefix)
    if not db_prefix_path.is_absolute():
        db_prefix_path = (sample_path / db_prefix_path).resolve()
    else:
        db_prefix_path = db_prefix_path.resolve()

    pin = db_prefix_path.parent / f"{db_prefix_path.name}.pin"
    psq = db_prefix_path.parent / f"{db_prefix_path.name}.psq"
    phr = db_prefix_path.parent / f"{db_prefix_path.name}.phr"
    idx_files = [p for p in (pin, psq, phr) if p.exists()]
    if len(idx_files) == 0:
        raise FileNotFoundError(f"BLAST DB index files not found for prefix: {db_prefix_path}")

    real_db_dir = idx_files[0].resolve().parent
    real_db_prefix = real_db_dir / db_prefix_path.name

    extra_mounts: List[Path] = []
    try:
        real_db_dir.relative_to(sample_path)
    except ValueError:
        extra_mounts.append(real_db_dir)

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

    # 1) 只有跑完 workspace init，Philosopher 才不会报 workspace not detected
    # 我们仅让它通过 init 和 annotate，再跑 iprophet 来合并，彻底抛弃 filter，防止崩溃
    _run(
        [
            "bash",
            "-lc",
            f"""
set -euo pipefail
cd {shlex.quote(str(work_dir))}
philosopher workspace --init
philosopher database --annotate {shlex.quote(str(combined_target_decoy_fasta))}
philosopher iprophet --threads {int(threads)} {shlex.quote(str(pep1))} {shlex.quote(str(pep2))} {shlex.quote(str(pep3))}
""",
        ],
        "philosopher-iprophet",
    )

    # 2) 用 Python 工具解析 iProphet 的结果文件，不惧怕任何截断或越界 Bug
    merged_xml = work_dir / "interact.iproph.pep.xml"
    _run(
        [
            "bash",
            "-lc",
            f"""
set -euo pipefail
cd {shlex.quote(str(work_dir))}
python /opt/tpp/bin/pepxml2csv.py {shlex.quote(str(merged_xml))}
""",
        ],
        "pepxml2csv",
    )

    merged_csv = work_dir / "interact.iproph.pep.csv"
    if not merged_csv.exists() and not dry_run:
        raise FileNotFoundError(f"Missing pepxml2csv output: {merged_csv}")

    # 3) 执行我们刚才写好的纯 Python 原生 FDR 过滤
    if not dry_run:
        fdr_csv = _filter_by_fdr(merged_csv, fdr)

        df_fdr = pd.read_csv(fdr_csv, sep="\t")
        pep_col = "peptide" if "peptide" in df_fdr.columns else ("Peptide" if "Peptide" in df_fdr.columns else None)
        if pep_col and not df_fdr.empty:
            peptides = sorted(set(df_fdr[pep_col].dropna().astype(str).tolist()))
        else:
            peptides = []

        blast_fa = work_dir / "peptide_blastp.fasta"
        blast_out = work_dir / "peptide_blastp.outfmt6.tsv"

        if peptides:
            with blast_fa.open("w", encoding="utf-8") as bf:
                for p in peptides:
                    bf.write(f">{p}\n{p}\n")

            # 4) Blastp 去除经典组蛋白 (Canonical removal)
            _run(
                [
                    "bash",
                    "-lc",
                    f"""
set -euo pipefail
cd {shlex.quote(str(work_dir))}
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
            df_final = df_fdr[~df_fdr[pep_col].astype(str).isin(remove_set)].copy()
        else:
            df_final = pd.DataFrame(columns=[pep_col] if pep_col else ["peptide"])
            if not blast_fa.exists():
                blast_fa.write_text("")
            if not blast_out.exists():
                blast_out.write_text("")

        # 5) MixMHCpred (如果输入了 HLA Alleles)
        if hla_alleles and not df_final.empty:
            bind_fa = work_dir / "binding.fasta"
            bind_out = work_dir / "mixmhcpred.out.tsv"

            bind_peptides = sorted(set([p for p in df_final[pep_col].astype(str).tolist() if 0 < len(p) < 15]))
            lines = []
            for p in bind_peptides:
                lines.append(f">{p}\n{p}\n")
            _write_text(bind_fa, "".join(lines))

            _run(
                [
                    "bash",
                    "-lc",
                    f"""
set -euo pipefail
cd {shlex.quote(str(work_dir))}
MixMHCpred -i {shlex.quote(str(bind_fa))} -o {shlex.quote(str(bind_out))} -a {shlex.quote(str(hla_alleles))}
""",
                ],
                "mixmhcpred",
            )

            if _nonempty(bind_out):
                rdf = pd.read_csv(bind_out, sep="\t", comment="#")
                if {"Peptide", "%Rank_bestAllele", "BestAllele"}.issubset(set(rdf.columns)):
                    sub = rdf[["Peptide", "%Rank_bestAllele", "BestAllele"]].drop_duplicates()
                    out_merge = df_final.merge(sub, left_on=pep_col, right_on="Peptide", how="left")
                    out_merge.drop(columns=["Peptide"], inplace=True, errors="ignore")
                else:
                    out_merge = df_final
            else:
                out_merge = df_final

            out_merge.to_csv(final_table, sep="\t", index=False)
        else:
            df_final.to_csv(final_table, sep="\t", index=False)