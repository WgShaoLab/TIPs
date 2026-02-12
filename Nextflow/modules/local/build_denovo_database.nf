nextflow.enable.dsl=2

process BUILD_DENOVO_DATABASE {
  tag "${sample_id}"

  input:
  tuple val(sample_id), path(blastp_dir), path(peptides_fasta)
  path te_npy
  path te_class_dic

  output:
  tuple val(sample_id), path("Denovo_TE_SoftMerged_InstaNovo.fasta"), emit: te_fasta
  tuple val(sample_id), path("Denovo_result_merged_tag_1.fasta"), emit: tagged_fasta

  script:
  """
  python << 'PYTHON_EOF'
import sys
sys.path.append('${params.code_base}')
import pandas as pd
import numpy as np
import os
from TE_immunopeptide.TE_searchPipline.Class.fasta_class import Fasta

peptide_soft_dic = {}
with open('${peptides_fasta}', 'r') as f:
    for line in f:
        if line.startswith('>'):
            seq = line[1:].strip()
            peptide_soft_dic[seq] = []

TE_pro_dic = np.load('${te_npy}', allow_pickle=True)[()]
TE_class_dic = np.load('${te_class_dic}', allow_pickle=True)[()]

blastp_dir = [d for d in os.listdir('.') if d.endswith('_Child_fasta')][0]
blastp_file = os.path.join(blastp_dir, 'merged_blastp_TE.csv')


df = pd.read_csv(blastp_file, sep='\\t')


if 'I/L_identity' in df.columns:
    df = df[df['I/L_identity'] == True]


df = df.drop_duplicates(subset=['qaccver', 'saccver'])



df['TE_name'] = df['saccver'].apply(lambda x: x.split('::')[0])
df['TE_class'] = df['TE_name'].apply(lambda x: TE_class_dic.get(x, 'Unknown'))

df_filter = df[['saccver', 'TE_class', 'TE_name']].drop_duplicates()

with open('Denovo_TE_SoftMerged_InstaNovo.fasta', 'w') as f:
    for _, row in df_filter.iterrows():
        te_header = row['saccver']
        seq = TE_pro_dic.get(te_header, '')
        if seq:
            te_class = row['TE_class']
            header = f"Denovo|mixed:::{row['TE_name']}:::{te_class}:::{te_header}"
            f.write(f'>{header}\\n{seq}\\n')

fasta = Fasta('Denovo_TE_SoftMerged_InstaNovo.fasta', uniform_PE_tag='1')
fasta.output('Denovo_result_merged_tag_1.fasta')
PYTHON_EOF
  """
}