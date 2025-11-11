import argparse
import os
import glob

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
from .utils import run_cmd, timing


# ----------------------------- helpers -----------------------------

def _find_files(patterns):
    """
    Collect files matching a list of glob patterns.
    Returns a de-duplicated, name-sorted list.
    """
    files = []
    for pat in patterns:
        files.extend(glob.glob(pat, recursive=True))
    # de-duplicate and keep order by using dict then back to list
    return sorted(list(dict.fromkeys(files)))


def _ensure_dir(path):
    """Create directory if not exists and return the path."""
    os.makedirs(path, exist_ok=True)
    return path


# ----------------------------- stages ------------------------------

@timing
def stage_denovo(cfg):
    """
    Run de novo engines, filter confident tags, and write the merged FASTA:
      {sample_path}/Denovo/Denovo_TE_SoftMerged_InstaNovo.fasta
    Header is the sequence itself. PE tags are unified later during DB building.
    """
    sample_path = cfg.sample_path
    gpu_id = cfg.gpu_id
    dcfg = cfg.raw["denovo"]

    filtered_tables = []

    # Casanovo
    if dcfg["casanovo"]["enable"]:
        eng = CasanovoEngine(
            sample_path=sample_path,
            gpu_id=gpu_id,
            binary=dcfg["casanovo"]["binary"],
            model_ckpt=dcfg["casanovo"]["model_ckpt"],
            config_yaml=dcfg["casanovo"]["config_yaml"],
        )
        eng.run()
        filtered_tables.append(
            eng.filter_subseq(dcfg["min_score"], dcfg["min_length"])
        )

    # PepNet
    if dcfg["pepnet"]["enable"]:
        eng = PepNetEngine(
            sample_path=sample_path,
            gpu_id=gpu_id,
            python_bin=dcfg["pepnet"]["python_bin"],
            script=dcfg["pepnet"]["script"],
            model_h5=dcfg["pepnet"]["model_h5"],
        )
        eng.run()
        filtered_tables.append(
            eng.filter_subseq(dcfg["min_score"], dcfg["min_length"])
        )

    # InstaNovo
    if dcfg["instanovo"]["enable"]:
        eng = InstaNovoEngine(
            sample_path=sample_path,
            gpu_id=gpu_id,
            python_bin=dcfg["instanovo"]["python_bin"],
            module_call=dcfg["instanovo"]["module_call"],
            model_ckpt=dcfg["instanovo"]["model_ckpt"],
            extra_env=dcfg["instanovo"].get("extra_env"),
        )
        eng.run()
        filtered_tables.append(
            eng.filter_subseq(dcfg["min_score"], dcfg["min_length"])
        )

    out_fa = os.path.join(sample_path, "Denovo", "Denovo_TE_SoftMerged_InstaNovo.fasta")
    write_denovo_merged_fasta(sample_path, filtered_tables, out_fa)


@timing
def stage_dbsearch(cfg):
    """
    Build search DBs and run DB search engines.
    The DB building stage already integrates:
      de novo tag → TE reference (blastp-short) → sample-specific TE sub-DB selection.
    """
    dbp = DBSearchPipeline(cfg)
    dbp.build_search_database()
    dbp.run_comet()
    dbp.run_msfragger()
    dbp.run_msgfplus()


@timing
def stage_postprocess(cfg):
    """
    Full post-processing workflow:
      1) Collect pepXML files → run RefreshParser on each
      2) Merge with xinteract and run iProphet (decoy prefix = rev_)
      3) Export table with pepxml2csv.py
      4) filter_by_fdr (charge-state aware FDR)
      5) blastp_remove_canonical (homology filtering against human proteome)
      6) predict_binding_affinity_mixmhcpred (HLA binding prediction)
    """
    sample_path = cfg.sample_path
    pcfg = cfg.raw["postprocessing"]
    tpp_cfg = pcfg["tpp"]

    # --- 1) Collect pepXML files ---
    # Accept both .pep.xml and .pepXML; search recursively under sample_path.
    pepxml_list = _find_files([
        os.path.join(sample_path, "**", "*.pep.xml"),
        os.path.join(sample_path, "**", "*.pepXML"),
    ])
    if len(pepxml_list) == 0:
        raise RuntimeError(
            "[[EN REQUIRED]] pepXML . pepXML;"
            "[[EN REQUIRED]] MSGF+ [[EN REQUIRED]] .mzid, pepXML ."
        )

    # Output locations
    iprophet_dir = _ensure_dir(os.path.join(sample_path, "DB_search_iProphet", "IPROPHET"))
    iprophet_prefix = os.path.join(iprophet_dir, "iprophet")
    iprophet_pepxml = f"{iprophet_prefix}.pep.xml"
    iprophet_csv = f"{iprophet_prefix}.tsv"   # pepxml2csv.py typically outputs TSV

    # --- 2) Run RefreshParser per pepXML ---
    refreshparser = tpp_cfg["refreshparser"]
    for pepxml in pepxml_list:
        cmd = f'"{refreshparser}" "{pepxml}"'
        run_cmd(cmd, env=os.environ)

    # --- 3) xinteract (merge + iProphet) ---
    # Use decoy prefix `rev_` (our decoy generator uses rev_)
    xinteract = tpp_cfg["xinteract"]
    # Common flags: -drev_ to specify decoy, -Op for default PP/iProphet behavior, -N for output prefix
    x_cmd = f'"{xinteract}" -drev_ -Op -N "{iprophet_prefix}" ' + " ".join(f'"{p}"' for p in pepxml_list)
    run_cmd(x_cmd, env=os.environ)

    # Some TPP builds prefer another RefreshParser pass on the merged pepXML
    if os.path.exists(iprophet_pepxml):
        run_cmd(f'"{refreshparser}" "{iprophet_pepxml}"', env=os.environ)
    else:
        # Fallback: locate the generated pepXML if the filename differs (e.g., *.interact.pep.xml)
        maybe = _find_files([os.path.join(iprophet_dir, "*.pep.xml")])
        if not maybe:
            raise RuntimeError("xinteract did not produce a merged pepXML; cannot continue.")
        iprophet_pepxml = maybe[0]
        run_cmd(f'"{refreshparser}" "{iprophet_pepxml}"', env=os.environ)

    # --- 4) pepxml2csv (export table) ---
    pepxml2csv = tpp_cfg["pepxml2csv"]
    # Try executing the script directly; fall back to `python script.py` if needed.
    try:
        run_cmd(f'"{pepxml2csv}" "{iprophet_pepxml}" > "{iprophet_csv}"', env=os.environ)
    except Exception:
        run_cmd(f'python "{pepxml2csv}" "{iprophet_pepxml}" > "{iprophet_csv}"', env=os.environ)

    # --- 5) FDR filtering ---
    fdr_thr = float(pcfg.get("fdr_threshold", 0.03))
    fdr_csv = filter_by_fdr(iprophet_csv, fdr_threshold=fdr_thr)

    # --- 6) Homology filtering against human proteome (blastp-short) ---
    # Only run when postprocessing.blastp.enable is True
    bcfg = pcfg.get("blastp", {})
    if bcfg.get("enable", True):
        blast_bin = bcfg["binary"]
        blast_db = bcfg["human_proteome_blastdb"]
        threads = int(bcfg.get("threads_per_job", 4))
        te_only_csv = blastp_remove_canonical(
            pep_csv=fdr_csv,
            blastp_bin=blast_bin,
            blast_db=blast_db,
            threads=threads,
        )
    else:
        te_only_csv = fdr_csv

    # --- 7) HLA binding prediction (MixMHCpred) ---
    hcfg = cfg.raw.get("hla_binding", {})
    if hcfg.get("enable", True):
        mixmhcpred_bin = hcfg["mixmhcpred_bin"]
        alleles = hcfg["hla_alleles"]
        pep_len_min = int(hcfg.get("peptide_length_min", 8))
        pep_len_max = int(hcfg.get("peptide_length_max", 14))
        final_csv = predict_binding_affinity_mixmhcpred(
            peptide_table=te_only_csv,
            mixmhcpred_bin=mixmhcpred_bin,
            hla_alleles=alleles,
            peptide_length_min=pep_len_min,
            peptide_length_max=pep_len_max,
        )
    else:
        final_csv = te_only_csv

    print("\n[stage_postprocess] Done:")
    print(f"  iProphet pepXML : {iprophet_pepxml}")
    print(f"  iProphet table   : {iprophet_csv}")
    print(f"  FDR-filtered    : {fdr_csv}")
    if te_only_csv != fdr_csv:
        print(f"  After homology removal    : {te_only_csv}")
    if final_csv != te_only_csv:
        print(f"  Binding prediction results    : {final_csv}")


# ----------------------------- CLI --------------------------------


def main():
    parser = argparse.ArgumentParser(description="TIPs pipeline CLI")
    parser.add_argument("--config", required=True, help="Path to YAML config")
    parser.add_argument(
        "--stage",
        required=True,
        choices=("denovo", "dbsearch", "postprocess", "all"),
        help="Which stage to run",
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
