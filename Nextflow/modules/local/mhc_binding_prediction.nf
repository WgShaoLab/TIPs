process MHC_BINDING_PREDICTION {
  tag "${sample_id}"
  container params.container_dbsearch
  publishDir "${params.outdir}/${sample_id}", mode: 'copy'

  input:
  tuple val(sample_id), path(filtered_peptides)

  output:
  tuple val(sample_id), path("final_TE_results.tsv"), emit: final_results
  tuple val(sample_id), path("binding.fasta"), emit: binding_fasta, optional: true
  tuple val(sample_id), path("binding.result"), emit: binding_result, optional: true

  script:
  """
  set -euo pipefail

  if [ ! -s "${filtered_peptides}" ]; then
    echo "Empty input file, copying as-is"
    cp "${filtered_peptides}" final_TE_results.tsv
    exit 0
  fi

  python3 << 'PYEOF'
import pandas as pd
import subprocess
import sys
import os

# Read filtered peptides
df = pd.read_csv("${filtered_peptides}", sep='\\t')

# Find peptide column
peptide_col = None
for c in df.columns:
    if c.lower() == 'peptide':
        peptide_col = c
        break

if peptide_col is None:
    print("No Peptide column found, copying input as-is")
    df.to_csv("final_TE_results.tsv", sep='\\t', index=False)
    sys.exit(0)

mhc_alleles = "${params.mhc_alleles}".strip()

# Match CLI behaviour: only run MixMHCpred if HLA alleles are explicitly provided
if not mhc_alleles:
    print("No HLA alleles specified, writing final results without MHC predictions")
    df.to_csv("final_TE_results.tsv", sep='\\t', index=False)
    sys.exit(0)

# Filter peptides for MHC-I length (8-14aa; legacy behavior: len < 15)
peptides = sorted(set([
    p for p in df[peptide_col].dropna().astype(str).tolist()
    if 0 < len(p) < 15
]))

print(f"Valid peptides for MHC prediction (len 8-14): {len(peptides)}")

if not peptides:
    print("No valid peptides for MHC prediction, copying input as-is")
    df.to_csv("final_TE_results.tsv", sep='\\t', index=False)
    sys.exit(0)

# Write binding FASTA
with open("binding.fasta", 'w') as f:
    for p in peptides:
        f.write(f">{p}\\n{p}\\n")

# Run MixMHCpred
print(f"Running MixMHCpred with alleles: {mhc_alleles}")

result = subprocess.run(
    ["MixMHCpred", "-i", "binding.fasta", "-o", "binding.result", "-a", mhc_alleles],
    capture_output=True, text=True
)
if result.returncode != 0:
    print(f"MixMHCpred failed (exit {result.returncode}): {result.stderr}")
    # Fallback: copy input as-is
    df.to_csv("final_TE_results.tsv", sep='\\t', index=False)
    sys.exit(0)

if not os.path.exists("binding.result") or os.path.getsize("binding.result") == 0:
    print("MixMHCpred produced no output, copying input as-is")
    df.to_csv("final_TE_results.tsv", sep='\\t', index=False)
    sys.exit(0)

# Parse MixMHCpred output
rdf = pd.read_csv("binding.result", sep='\\t', comment='#')

# Robust column detection for different MixMHCpred output formats
if {"Peptide", "%Rank_bestAllele", "BestAllele"}.issubset(set(rdf.columns)):
    sub = rdf[["Peptide", "%Rank_bestAllele", "BestAllele"]].drop_duplicates()
    # Merge using pandas (exact matching, no regex issues); preserve all input columns
    out = df.merge(sub, left_on=peptide_col, right_on="Peptide", how="left")
    out.drop(columns=["Peptide"], inplace=True, errors="ignore")
else:
    print("Warning: MixMHCpred output missing expected columns, copying input as-is")
    out = df

out.to_csv("final_TE_results.tsv", sep='\\t', index=False)

# Summary
total = len(out)
with_pred = out["%Rank_bestAllele"].notna().sum() if "%Rank_bestAllele" in out.columns else 0
print(f"Final results: {with_pred}/{total} peptides with MHC predictions")
PYEOF
  """
}
