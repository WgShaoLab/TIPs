nextflow.enable.dsl=2

process IPROPHET_INTEGRATION_TE {
  tag "${sample_id}_${type}"

  input:
  tuple val(sample_id), val(type), path(comet_pepxml), path(msfragger_pepxml), path(msgf_pepxml), path(database)

  output:
  tuple val(sample_id), val(type), path("interact.iproph.pep.xml"), emit: iprophet_pepxml
  tuple val(sample_id), path("peptide.tsv"), emit: peptide_tsv

  script:
  """
  set -euo pipefail

  philosopher workspace --init
  philosopher database --annotate ${database}

  /opt/tpp/bin/RefreshParser ${comet_pepxml} ${database}
  /opt/tpp/bin/RefreshParser ${msfragger_pepxml} ${database}
  /opt/tpp/bin/RefreshParser ${msgf_pepxml} ${database}

  philosopher iprophet --threads ${params.threads} \\
    ${comet_pepxml} ${msfragger_pepxml} ${msgf_pepxml}

  python3 << 'PYEOF'
import xml.etree.ElementTree as ET
import pandas as pd
from collections import defaultdict
import sys, os

# -------------------------------------------------------------------
# 1. Parse iProphet pepXML -> PSM-level DataFrame
#    (mirrors CLI integrate_te._parse_iproph_pepxml_to_psm_df)
# -------------------------------------------------------------------
def _is_decoy_protein(protein_id):
    pid = (protein_id or "").strip()
    return pid.startswith("rev_") or pid.startswith("DECOY_")

def parse_iproph_pepxml(pepxml_path):
    rows = []
    context = ET.iterparse(pepxml_path, events=("start", "end"))
    _, root = next(context)

    current_spec = None
    current_charge = None

    for event, elem in context:
        tag = elem.tag.split("}")[-1] if "}" in elem.tag else elem.tag

        if event == "start" and tag == "spectrum_query":
            current_spec = elem.attrib.get("spectrum") or elem.attrib.get("spectrumNativeID")
            current_charge = elem.attrib.get("assumed_charge")

        if event == "end" and tag == "search_hit":
            peptide = elem.attrib.get("peptide")
            protein = elem.attrib.get("protein")

            proteins = []
            if protein:
                proteins.append(protein)

            for ap in elem.findall(".//{*}alternative_protein"):
                ap_prot = ap.attrib.get("protein")
                if ap_prot:
                    proteins.append(ap_prot)

            iprob = None
            for ip in elem.findall(".//{*}interprophet_result"):
                p = ip.attrib.get("probability")
                if p is not None:
                    try:
                        iprob = float(p)
                    except Exception:
                        iprob = None

            if iprob is None:
                for pp in elem.findall(".//{*}peptideprophet_result"):
                    p = pp.attrib.get("probability")
                    if p is not None:
                        try:
                            iprob = float(p)
                        except Exception:
                            iprob = None

            if peptide and (iprob is not None) and current_spec:
                uniq_proteins = sorted(set(proteins)) if proteins else []
                is_decoy = True
                if uniq_proteins:
                    is_decoy = all(_is_decoy_protein(x) for x in uniq_proteins)

                rows.append({
                    "spectrum": current_spec,
                    "assumed_charge": int(current_charge) if current_charge else None,
                    "peptide": peptide,
                    "proteins": ";".join(uniq_proteins),
                    "is_decoy": bool(is_decoy),
                    "iprob": float(iprob),
                })

            elem.clear()

        if event == "end" and tag == "spectrum_query":
            elem.clear()

    root.clear()
    return pd.DataFrame(rows)


# -------------------------------------------------------------------
# 2. Compute q-values (monotone-transformed cumulative FDR)
#    (mirrors CLI integrate_te._compute_qvalues_from_sorted)
# -------------------------------------------------------------------
def compute_qvalues_from_sorted(is_decoy):
    fdrs = []
    d = 0
    t = 0
    for flag in is_decoy:
        if flag:
            d += 1
        else:
            t += 1
        if t == 0:
            fdrs.append(1.0)
        else:
            fdrs.append(d / t)

    qvals = [0.0] * len(fdrs)
    running_min = 1.0
    for i in range(len(fdrs) - 1, -1, -1):
        running_min = min(running_min, fdrs[i])
        qvals[i] = running_min
    return qvals


# -------------------------------------------------------------------
# 3. sp| filter: remove peptides mapped to canonical UniProt proteins
#    (mirrors CLI integrate_te._filter_sp_mapped_peptides)
# -------------------------------------------------------------------
def filter_sp_mapped_peptides(df, col="Mapped Proteins"):
    mask = df[col].astype(str).str.contains(r"sp\\|", regex=True)
    return df[~mask].copy()


# ===================================================================
# MAIN
# ===================================================================
pepxml = "interact.iproph.pep.xml"
fdr_threshold = float("${params.fdr_threshold}")

# Parse
psm_all = parse_iproph_pepxml(pepxml)
if psm_all.empty:
    print("ERROR: No PSMs parsed from iProphet pepXML", file=sys.stderr)
    sys.exit(1)

# Sort by iprob descending
psm_all = psm_all.sort_values(["iprob"], ascending=False).reset_index(drop=True)

# --- PSM-level FDR ---
psm_all["q_value"] = compute_qvalues_from_sorted(psm_all["is_decoy"].astype(bool).tolist())
psm_f = psm_all[psm_all["q_value"] <= fdr_threshold].copy()
psm_f.to_csv("psm.tsv", sep="\\t", index=False)

# --- Peptide-level FDR ---
# Collapse by peptide, use best iprob
pep_best = (
    psm_all.sort_values(["iprob"], ascending=False)
    .groupby("peptide", as_index=False)
    .first()
)
pep_best = pep_best.sort_values(["iprob"], ascending=False).reset_index(drop=True)
pep_best["q_value"] = compute_qvalues_from_sorted(pep_best["is_decoy"].astype(bool).tolist())
pep_f = pep_best[pep_best["q_value"] <= fdr_threshold].copy()

# Rename columns to match philosopher-like output
pep_f.rename(columns={"proteins": "Mapped Proteins", "iprob": "Probability"}, inplace=True)

# --- sp| filter: remove canonical UniProt-mapped peptides ---
pep_f = filter_sp_mapped_peptides(pep_f)

# Write all columns to match CLI behaviour (spectrum, assumed_charge, is_decoy, etc.)
pep_f.to_csv("peptide.tsv", sep="\\t", index=False)

n_psm = len(psm_f)
n_pep = len(pep_f)
print(f"FDR-filtered PSMs: {n_psm}")
print(f"FDR-filtered peptides (after sp| filter): {n_pep}")
PYEOF
  """
}