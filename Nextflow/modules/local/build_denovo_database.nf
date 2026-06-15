nextflow.enable.dsl=2

process BUILD_DENOVO_DATABASE {
  tag "${sample_id}"

  input:
  tuple val(sample_id), path(blastp_dir), path(peptides_fasta), path(peptide_soft_json)
  path te_npy
  path te_class_dic

  output:
  tuple val(sample_id), path("Denovo_TE_SoftMerged_Merged.fasta"), emit: te_fasta

  script:
  """
  python << 'PYTHON_EOF'
import sys
sys.path.append('${params.code_base}')
import json
import pandas as pd
import numpy as np
import os

with open('${peptide_soft_json}', 'r') as f:
    peptide_soft_dic = json.load(f)

te_pro_dic = np.load('${te_npy}', allow_pickle=True).item()
te_class_dic = np.load('${te_class_dic}', allow_pickle=True).item()

blastp_dir = [d for d in os.listdir('.') if d.endswith('_Child_fasta')][0]
blastp_file = os.path.join(blastp_dir, 'merged_blastp_TE.csv')

df = pd.read_csv(blastp_file, sep='\\t')

if 'I/L_identity' in df.columns:
    df = df[df['I/L_identity'] == True]

# Build soft column from peptide_soft_dic (match CLI denovo_te.py)
df['soft'] = df['qaccver'].astype(str).apply(lambda x: peptide_soft_dic.get(x, []))

# TE_name extraction matching CLI _te_name_from_saccver
def _te_name_from_saccver(x):
    s = str(x)
    if "Class=DNA;" not in s:
        return s.split("::")[0]
    return s.split(";")[0].split("=")[1]

df['TE_name'] = df['saccver'].astype(str).apply(_te_name_from_saccver)
df['TE_class'] = df['TE_name'].astype(str).apply(lambda x: te_class_dic[x] if x in te_class_dic else "DNA")

df = df.drop_duplicates(subset=['qaccver', 'saccver'])

# Group and aggregate soft tags (match CLI denovo_te.py:698-702)
out_df = df[['saccver', 'TE_class', 'soft', 'TE_name']].copy()
out_df = (
    out_df.groupby(['saccver', 'TE_class', 'TE_name'], as_index=False)
    .agg({'soft': lambda xs: sum(xs, [])})
)
out_df['soft'] = out_df['soft'].apply(lambda xs: "_".join(sorted(set(xs))))

with open('Denovo_TE_SoftMerged_Merged.fasta', 'w') as f:
    for _, row in out_df.iterrows():
        soft = row['soft']
        te_header = row['saccver']
        te_name = row['TE_name']
        te_class = row['TE_class']
        if te_header not in te_pro_dic:
            continue
        seq = te_pro_dic[te_header]
        header = f"Denovo|{soft}:::{te_name}:::{te_class}:::{te_header}"
        f.write(f'>{header}\\n{seq}\\n')

PYTHON_EOF
  """
}