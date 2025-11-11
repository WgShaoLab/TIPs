import os
import re
import ast
import time
import numpy as np
import pandas as pd
import subprocess
from typing import List, Dict
from .utils import run_cmd, wait_while_running


def filter_by_fdr(
    peptide_csv: str,
    fdr_threshold: float,
    out_path: str | None = None
) -> str:
    """
    Compute charge-specific running FDR (decoy/target),
    keep only PSMs passing fdr_threshold.
    Assumes the CSV is PeptideProphet / iProphet-like export
    with columns including:
        - analysis_result (list-like string containing peptideprophet_result.probability)
        - protein
        - assumed_charge
    Decoys are prefixed with "rev_".
    """
    if out_path is None:
        out_path = peptide_csv.replace(".pep.csv", "_fdr.csv")

    df_merged = pd.DataFrame()
    try:
        df_all = pd.read_csv(peptide_csv, sep="\t")
    except Exception:
        print("[filter_by_fdr] CSV appears empty or missing.")
        df_merged.to_csv(out_path, index=False, sep="\t")
        return out_path

    if len(df_all) < 20:
        df_merged.to_csv(out_path, index=False, sep="\t")
        return out_path

    # primary protein (strip multi-assignments after '#')
    df_all["protein_primary"] = df_all["protein"].apply(lambda x: x.split("#")[0])

    # parse peptideProphet probability
    df_all["peptideProphet"] = df_all["analysis_result"].apply(ast.literal_eval)
    df_all["peptideProphet_probability"] = df_all["peptideProphet"].apply(
        lambda x: float(x[0]["peptideprophet_result"]["probability"])
    )

    # mark decoy
    df_all["decoy"] = df_all["protein_primary"].apply(
        lambda x: 1 if "rev_" in x else 0
    )

    # sort by prob desc
    df_all = df_all.sort_values(
        by=["peptideProphet_probability"],
        ascending=False
    ).reset_index(drop=True)

    # process each charge state separately
    for charge, df_c in df_all.groupby("assumed_charge"):
        df_c = df_all[df_all["assumed_charge"] == charge].copy()

        # cumulative counts
        df_c["target_cumsum"] = (df_c["decoy"] == 0).cumsum()
        df_c["decoy_cumsum"] = (df_c["decoy"] == 1).cumsum()

        # naive FDR ~ decoy/target
        df_c["fdr"] = df_c["decoy_cumsum"] / df_c["target_cumsum"]

        reversed_df = df_c[::-1]
        valid_idx = reversed_df[reversed_df["fdr"] < fdr_threshold].index
        if len(valid_idx) == 0:
            print(f"[filter_by_fdr] charge {charge}: no PSM with FDR < {fdr_threshold}")
            continue
        keep_until = valid_idx[0]
        filtered_df = df_c.iloc[: keep_until + 1]
        filtered_df = filtered_df[filtered_df["decoy"] == 0]  # keep only targets

        print(f"[filter_by_fdr] charge {charge}: {len(filtered_df)} kept (FDR<{fdr_threshold})")
        df_merged = pd.concat([df_merged, filtered_df]) if not df_merged.empty else filtered_df

    if df_merged.empty:
        df_merged.to_csv(out_path, index=False, sep="\t")
    else:
        df_merged = df_merged.sort_values(
            by=["peptideProphet_probability"],
            ascending=False
        )
        df_merged.to_csv(out_path, index=False, sep="\t")

    return out_path


def blastp_remove_canonical(
    pep_csv: str,
    blastp_bin: str,
    blast_db: str,
    threads: int = 20,
    poll_interval: int = 5,
) -> str:
    """
    Remove canonical peptides by blasting peptides against a
    reference (canonical) proteome DB.
    Logic:
      - Create a special FASTA containing each peptide as both header+seq.
      - Run blastp-short.
      - Keep only peptides that are NOT perfect full-length canonical matches.
    """
    df = pd.read_csv(pep_csv, sep="\t")
    if "peptide" in df.columns:
        pep_col = "peptide"
    elif "Peptide" in df.columns:
        pep_col = "Peptide"
    else:
        raise ValueError("Input CSV must contain 'peptide' or 'Peptide' column")

    unique_peps = sorted(set(df[pep_col].tolist()))

    fasta_for_blast = pep_csv.replace(".csv", "_blastp.fasta").replace(".tsv", "_blastp.fasta")
    with open(fasta_for_blast, "w") as fh:
        for seq in unique_peps:
            fh.write(f">{seq}\n{seq}\n")

    blast_out = pep_csv.replace(".csv", "_blastpResult.txt").replace(".tsv", "_blastpResult.txt")
    cmd = (
        f"{blastp_bin} -task blastp-short "
        f"-query {fasta_for_blast} "
        f"-db {blast_db} "
        f"-out {blast_out} "
        f"-outfmt 6 -evalue 20000 -num_threads {threads} &"
    )
    run_cmd(cmd)

    # wait for blastp to finish
    while True:
        try:
            n = int(
                subprocess.check_output(
                    "ps aux | grep -v grep | grep -c 'blastp'",
                    shell=True,
                ).strip()
            )
        except Exception:
            n = 0
        if n <= 0:
            break
        time.sleep(poll_interval)

    # parse blast results
    cols = [
        "Identifier",
        "Protein",
        "Score",
        "Length",
        "Start",
        "End",
        "MatchStart",
        "MatchEnd",
        "DbStart",
        "DbEnd",
        "Value1",
        "Value2",
    ]
    bdf = pd.read_csv(blast_out, sep="\t", names=cols)
    bdf["pep_len"] = bdf["Identifier"].apply(len)

    # keep only perfect full-length 100% identity alignments
    full_hits = bdf[
        (bdf["MatchStart"] == 1) &
        (bdf["MatchEnd"] == bdf["pep_len"]) &
        (bdf["Score"] == 100)
    ]
    remove_list = full_hits["Identifier"].drop_duplicates().tolist()

    df_filtered = df[~df[pep_col].isin(remove_list)]
    filtered_out = pep_csv.replace(".csv", "_blastp.filter.txt").replace(".tsv", "_blastp.filter.txt")
    df_filtered.to_csv(filtered_out, sep="\t", index=False)

    print("[blastp_remove_canonical] canonical-like peptides removed")
    return filtered_out


def predict_binding_affinity_mixmhcpred(
    peptide_table: str,
    mixmhcpred_bin: str,
    hla_alleles: str,
    peptide_length_min: int = 8,
    peptide_length_max: int = 14,
    out_csv: str | None = None,
) -> str:
    """
    Predict HLA binding affinity using MixMHCpred.
    Only peptides in [peptide_length_min, peptide_length_max] are tested.

    peptide_table:
      TSV with column "Peptide" (after iProphet/blastp filtering).
    """
    df = pd.read_csv(peptide_table, sep="\t")
    if len(df) == 0:
        print("[predict_binding_affinity_mixmhcpred] no peptides, skipping.")
        return peptide_table

    binding_fasta = os.path.join(os.path.dirname(peptide_table), "binding.fasta")
    if os.path.exists(binding_fasta):
        os.remove(binding_fasta)

    with open(binding_fasta, "w") as fh:
        for pep in df["Peptide"]:
            if peptide_length_min < len(pep) < peptide_length_max:
                fh.write(f">{pep}\n{pep}\n")

    out_result = binding_fasta.replace(".fasta", ".result")
    run_cmd(
        f"{mixmhcpred_bin} "
        f"-i {binding_fasta} "
        f"-o {out_result} "
        f"-a {hla_alleles}"
    )

    pred_df = pd.read_csv(out_result, sep="\t", comment="#")
    pred_df = pred_df[["Peptide", "%Rank_bestAllele", "BestAllele"]]

    rank_map = dict(zip(pred_df["Peptide"], pred_df["%Rank_bestAllele"]))
    allele_map = dict(zip(pred_df["Peptide"], pred_df["BestAllele"]))

    df["%Rank_bestAllele"] = df["Peptide"].map(rank_map)
    df["BestAllele"] = df["Peptide"].map(allele_map).fillna("")

    if out_csv is None:
        out_csv = os.path.join(os.path.dirname(peptide_table), "final_TE_results.csv")
    df.to_csv(out_csv, sep="\t", index=False)
    print(f"[predict_binding_affinity_mixmhcpred] wrote {out_csv}")
    return out_csv
