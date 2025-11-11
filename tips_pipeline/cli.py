import argparse
import os
import glob
import subprocess
import sys
import shutil
from typing import List

import pandas as pd  # needed for dual-track splitting

from .config import load_config
from .pipelines_denovo import (
    CasanovoEngine,
    PepNetEngine,
    InstaNovoEngine,
    write_denovo_merged_fasta,
)
from .pipelines_dbsearch import DBSearchPipeline
from .ms_postprocess import (
    filter_by_fdr,
    blastp_remove_canonical,
    predict_binding_affinity_mixmhcpred,
)


# -------------------------- small helpers --------------------------

def run_cmd(cmd: str, env=None):
    """
    Run a shell command and stream output. Raise on non-zero return code.
    """
    print(f"[CMD] {cmd}")
    ret = subprocess.run(cmd, shell=True, env=env)
    if ret.returncode != 0:
        raise RuntimeError(f"Command failed with code {ret.returncode}: {cmd}")


def _ensure_dir(path: str) -> str:
    """
    Create directory if missing and return the path.
    """
    os.makedirs(path, exist_ok=True)
    return path


def _first_exists(paths: List[str]) -> str:
    """
    Return the first path that exists in the list, else empty string.
    """
    for p in paths:
        if p and os.path.exists(p):
            return p
    return ""


# -------------------------- stages --------------------------

def stage_denovo(cfg):
    """
    Stage 1: Run de novo engines based on config and write a merged de novo FASTA.

    This keeps the structure that other parts of the pipeline expect.
    """
    sample_path = cfg.sample_path
    sample_name = cfg.raw["sample"]["name"]

    denovo_cfg = cfg.raw.get("denovo", {})
    enable = bool(denovo_cfg.get("enable", False))
    if not enable:
        print("[stage_denovo] Skipped (disabled in config).")
        return

    engines = denovo_cfg.get("engines", [])
    out_dir = _ensure_dir(os.path.join(sample_path, "DeNovo"))
    fasta_out = os.path.join(out_dir, f"{sample_name}.denovo.fasta")

    # Collect per-engine peptide tables or FASTAs, then write merged FASTA.
    collected = []

    if "casanovo" in engines:
        print("[stage_denovo] Running Casanovo...")
        engine = CasanovoEngine(cfg)
        collected.append(engine.run_and_collect())

    if "pepnet" in engines:
        print("[stage_denovo] Running PepNet...")
        engine = PepNetEngine(cfg)
        collected.append(engine.run_and_collect())

    if "instanovo" in engines:
        print("[stage_denovo] Running InstaNovo...")
        engine = InstaNovoEngine(cfg)
        collected.append(engine.run_and_collect())

    print("[stage_denovo] Writing merged de novo FASTA...")
    write_denovo_merged_fasta(collected, fasta_out)
    print(f"[stage_denovo] Done. Merged de novo FASTA: {fasta_out}")


def stage_dbsearch(cfg):
    """
    Stage 2: Database search against (proteome +/- de novo) FASTA using external engines
    configured in DBSearchPipeline, producing pepXMLs for each search.
    """
    enable = bool(cfg.raw.get("dbsearch", {}).get("enable", True))
    if not enable:
        print("[stage_dbsearch] Skipped (disabled in config).")
        return

    pipe = DBSearchPipeline(cfg)
    pipe.run()
    print("[stage_dbsearch] Done.")


def stage_postprocess(cfg):
    """
    Stage 3: TPP-based iProphet merge; then either:
      - single_merged: one FDR + (optional) homology removal + (optional) HLA prediction
      - dual_tracks : split by source tags (SRC=TE / SRC=HUMAN) and perform independent QC.
                      In dual_tracks, both TE and Canonical tracks undergo independent FDR.
                      TE then optionally runs homology removal and HLA prediction; Canonical
                      can also run optional HLA prediction (controlled by hla_binding.enable).
    """
    sample_path = cfg.sample_path
    sample_name = cfg.raw["sample"]["name"]

    pcfg = cfg.raw.get("postprocessing", {})
    if not pcfg.get("enable", True):
        print("[stage_postprocess] Skipped (disabled in config).")
        return

    tpp_cfg = pcfg.get("tpp", {})
    ipcfg = pcfg.get("iprophet", {})

    # Output locations
    iprophet_dir = _ensure_dir(os.path.join(sample_path, "DB_search_iProphet", "IPROPHET"))
    iprophet_prefix = os.path.join(iprophet_dir, "iprophet")
    iprophet_pepxml = f"{iprophet_prefix}.pep.xml"
    iprophet_csv = f"{iprophet_prefix}.tsv"  # pepxml2csv typically outputs TSV

    # --- 1) collect pepXMLs (from DB search) ---
    # We try common folders; fallback to recursive search.
    pepxml_roots = [
        os.path.join(sample_path, "DB_search_iProphet"),
        os.path.join(sample_path, "DB_search"),
        os.path.join(sample_path, "Search"),
    ]
    pepxml_list = []
    for root in pepxml_roots:
        if os.path.isdir(root):
            pepxml_list.extend(glob.glob(os.path.join(root, "**", "*.pep.xml"), recursive=True))
            pepxml_list.extend(glob.glob(os.path.join(root, "**", "*.pepXML"), recursive=True))
    # de-duplicate
    pepxml_list = sorted(set(pepxml_list))

    if not pepxml_list:
        # final fallback to sample_path
        pepxml_list = sorted(set(
            glob.glob(os.path.join(sample_path, "**", "*.pep.xml"), recursive=True) +
            glob.glob(os.path.join(sample_path, "**", "*.pepXML"), recursive=True)
        ))

    if not pepxml_list:
        raise FileNotFoundError("[stage_postprocess] No pepXML files found for iProphet.")

    print(f"[stage_postprocess] Found {len(pepxml_list)} pepXML files.")

    # --- 2) RefreshParser per pepXML (add required attributes for TPP) ---
    refreshparser = tpp_cfg.get("refreshparser", "")
    if not refreshparser:
        raise ValueError("[stage_postprocess] Missing 'postprocessing.tpp.refreshparser' in config.")
    for pepxml in pepxml_list:
        run_cmd(f'"{refreshparser}" "{pepxml}"', env=os.environ)

    # --- 3) xinteract (merge + iProphet) ---
    # We use decoy prefix 'rev_' by default (consistent with typical reversed decoy generation).
    xinteract = tpp_cfg.get("xinteract", "")
    if not xinteract:
        raise ValueError("[stage_postprocess] Missing 'postprocessing.tpp.xinteract' in config.")

    # Common flags example: -drev_ -OIPROPHET/path/iprophet -Niprophet
    # -N sets output prefix; -O ensures output path; -d sets decoy prefix.
    decoy_prefix = pcfg.get("decoy_prefix", "rev_")
    interact_cmd = (
        f'"{xinteract}" -d{decoy_prefix} -O"{iprophet_dir}" -N"{os.path.basename(iprophet_prefix)}" '
        + " ".join([f'"{p}"' for p in pepxml_list])
    )
    run_cmd(interact_cmd, env=os.environ)

    # --- 4) pepxml2csv (export iProphet table) ---
    pepxml2csv = tpp_cfg.get("pepxml2csv", "")
    if not pepxml2csv:
        raise ValueError("[stage_postprocess] Missing 'postprocessing.tpp.pepxml2csv' in config.")
    # Try as executable; fallback to "python <script>" if needed.
    try:
        run_cmd(f'"{pepxml2csv}" "{iprophet_pepxml}" > "{iprophet_csv}"', env=os.environ)
    except Exception:
        run_cmd(f'{sys.executable} "{pepxml2csv}" "{iprophet_pepxml}" > "{iprophet_csv}"', env=os.environ)

    # --------------------- dual-track processing ---------------------
    mode = str(ipcfg.get("mode", "single_merged")).strip().lower()
    tracks = [t.strip().lower() for t in ipcfg.get("tracks", ["TE", "Canonical"])]

    if mode == "dual_tracks":
        print("[stage_postprocess] dual_tracks mode enabled. Splitting by source tags in FASTA headers...")
        # 2.1 Load the merged iProphet table once
        if not os.path.exists(iprophet_csv):
            raise FileNotFoundError(f"[stage_postprocess] iProphet table missing: {iprophet_csv}")
        df_all = pd.read_csv(iprophet_csv, sep="\t").fillna("")

        # 2.2 Build masks using source tags embedded in protein headers during database construction.
        #      TE entries should contain 'SRC=TE'
        #      Canonical (human proteome) entries should contain 'SRC=HUMAN'
        protein_col = "protein" if "protein" in df_all.columns else _first_exists(
            [c for c in df_all.columns if c.lower() in ("protein", "proteins", "accession", "accessions")]
        )
        if not protein_col:
            raise KeyError("[stage_postprocess] Could not locate a protein/accession column in iProphet table.")

        series = df_all[protein_col].astype(str)

        # Exclude contaminants explicitly if labeled
        contam_mask = series.str.contains("SRC=CONTAM", na=False)

        te_mask = series.str.contains("SRC=TE", na=False)
        human_mask = series.str.contains("SRC=HUMAN", na=False)

        te_df = df_all[te_mask & (~contam_mask)].copy()
        can_df = df_all[human_mask & (~contam_mask)].copy()

        # 2.3 Write per-track raw tables for traceability
        te_tsv = os.path.join(iprophet_dir, "iprophet_TE.tsv")
        can_tsv = os.path.join(iprophet_dir, "iprophet_Canonical.tsv")
        te_df.to_csv(te_tsv, sep="\t", index=False)
        can_df.to_csv(can_tsv, sep="\t", index=False)
        print(f"[stage_postprocess] Wrote track tables:\n  TE: {te_tsv}\n  Canonical: {can_tsv}")

        # 2.4 Independent FDR filtering for each track
        fdr_thr = float(pcfg.get("fdr_threshold", 0.03))
        fdr_te_csv = filter_by_fdr(te_tsv, fdr_threshold=fdr_thr)
        fdr_can_csv = filter_by_fdr(can_tsv, fdr_threshold=fdr_thr)
        print(f"[stage_postprocess] FDR threshold={fdr_thr:.3f}")
        print(f"[stage_postprocess] TE FDR table      : {fdr_te_csv}")
        print(f"[stage_postprocess] Canonical FDR table: {fdr_can_csv}")

        # 2.5 TE track: homology removal against human proteome (blastp-short), if enabled
        bcfg = pcfg.get("blastp", {})
        te_after_homology = fdr_te_csv
        if bcfg.get("enable", True):
            te_after_homology = blastp_remove_canonical(
                pep_csv=fdr_te_csv,
                blastp_bin=bcfg["binary"],
                blast_db=bcfg["human_proteome_blastdb"],
                threads=int(bcfg.get("threads_per_job", 4)),
            )
            print(f"[stage_postprocess] TE after homology removal: {te_after_homology}")

        # 2.6 Optional HLA binding prediction (both tracks if enabled)
        hcfg = cfg.raw.get("hla_binding", {})
        final_te_csv = te_after_homology
        if hcfg.get("enable", True):
            final_te_csv = predict_binding_affinity_mixmhcpred(
                peptide_table=te_after_homology,
                mixmhcpred_bin=hcfg["mixmhcpred_bin"],
                hla_alleles=hcfg["hla_alleles"],
                peptide_length_min=int(hcfg.get("peptide_length_min", 8)),
                peptide_length_max=int(hcfg.get("peptide_length_max", 14)),
            )
            print(f"[stage_postprocess] TE binding prediction results: {final_te_csv}")

        final_can_csv = fdr_can_csv
        if hcfg.get("enable", True):
            final_can_csv = predict_binding_affinity_mixmhcpred(
                peptide_table=fdr_can_csv,
                mixmhcpred_bin=hcfg["mixmhcpred_bin"],
                hla_alleles=hcfg["hla_alleles"],
                peptide_length_min=int(hcfg.get("peptide_length_min", 8)),
                peptide_length_max=int(hcfg.get("peptide_length_max", 14)),
            )
            print(f"[stage_postprocess] Canonical binding prediction results: {final_can_csv}")

        print("\n[stage_postprocess / dual_tracks] Done:")
        print(f"  iProphet pepXML               : {iprophet_pepxml}")
        print(f"  iProphet table (all)          : {iprophet_csv}")
        print(f"  Canonical final table         : {final_can_csv}")
        print(f"  TE final table                : {final_te_csv}")
        return

    # --------------------- original single-merged processing ---------------------
    # --- 5) FDR filtering (one table) ---
    fdr_thr = float(pcfg.get("fdr_threshold", 0.03))
    merged_fdr_csv = filter_by_fdr(iprophet_csv, fdr_threshold=fdr_thr)
    print(f"[stage_postprocess] single_merged FDR threshold={fdr_thr:.3f}")
    print(f"[stage_postprocess] FDR table: {merged_fdr_csv}")

    # --- 6) Optional homology filtering on merged table ---
    bcfg = pcfg.get("blastp", {})
    after_homology_csv = merged_fdr_csv
    if bcfg.get("enable", True):
        after_homology_csv = blastp_remove_canonical(
            pep_csv=merged_fdr_csv,
            blastp_bin=bcfg["binary"],
            blast_db=bcfg["human_proteome_blastdb"],
            threads=int(bcfg.get("threads_per_job", 4)),
        )
        print(f"[stage_postprocess] After homology removal: {after_homology_csv}")

    # --- 7) Optional HLA binding prediction on merged table ---
    hcfg = cfg.raw.get("hla_binding", {})
    final_csv = after_homology_csv
    if hcfg.get("enable", True):
        final_csv = predict_binding_affinity_mixmhcpred(
            peptide_table=after_homology_csv,
            mixmhcpred_bin=hcfg["mixmhcpred_bin"],
            hla_alleles=hcfg["hla_alleles"],
            peptide_length_min=int(hcfg.get("peptide_length_min", 8)),
            peptide_length_max=int(hcfg.get("peptide_length_max", 14)),
        )
        print(f"[stage_postprocess] Binding prediction results: {final_csv}")

    print("\n[stage_postprocess / single_merged] Done:")
    print(f"  iProphet pepXML  : {iprophet_pepxml}")
    print(f"  iProphet table   : {iprophet_csv}")
    print(f"  Final table      : {final_csv}")


# -------------------------- CLI --------------------------

def main():
    parser = argparse.ArgumentParser(
        description="TIPsSearch CLI: de novo, DB search, and post-processing (iProphet)."
    )
    parser.add_argument(
        "-c", "--config", required=True, help="Path to RUN_TIPs.yaml (or equivalent) config."
    )
    parser.add_argument(
        "-s", "--stage", choices=["denovo", "dbsearch", "postprocess", "all"], default="all",
        help="Which stage to run."
    )
    args = parser.parse_args()

    cfg = load_config(args.config)

    if args.stage in ("denovo", "all"):
        stage_denovo(cfg)
    if args.stage in ("dbsearch", "all"):
        stage_dbsearch(cfg)
    if args.stage in ("postprocess", "all"):
        stage_postprocess(cfg)


if __name__ == "__main__":
    main()
