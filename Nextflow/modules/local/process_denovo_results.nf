nextflow.enable.dsl=2

process PROCESS_DENOVO_RESULTS {
  tag "${sample_id}"

  input:
  tuple val(sample_id), path(casanovo_files), path(pepnet_files), path(instanovo_files)
  path te_npy
  path te_class_dic

  output:
  tuple val(sample_id), path("merged_raw_peptide.fasta"), emit: peptides

  script:
  """
  set -euo pipefail
  export CUDA_VISIBLE_DEVICES=""

  mkdir -p CasaNovo Pepnet Instanovo

  cp ${casanovo_files} CasaNovo/ 2>/dev/null || true
  cp ${pepnet_files} Pepnet/ 2>/dev/null || true
  cp ${instanovo_files} Instanovo/ 2>/dev/null || true

  python << 'PYTHON_EOF'
import sys
sys.path.append('${params.code_base}')

import os
os.environ['CUDA_VISIBLE_DEVICES'] = ''

import pandas as pd
import numpy as np
from io import StringIO
import re

from TE_immunopeptide.TE_searchPipline.SoftPipeline_class_iProphet import extract_subsequences_sub

def process_pepnet(path, min_score=0.75, min_length=8):
    file_list = [x for x in os.listdir(path) if x.endswith('.result')]
    df = pd.DataFrame()
    for file in file_list:
        fpath = os.path.join(path, file)
        if df.empty:
            df = pd.read_csv(fpath, sep='\\t')
        else:
            df = pd.concat([df, pd.read_csv(fpath, sep='\\t')])

    filter_df = df[['DENOVO', 'Positional Score']].copy()
    filter_df = filter_df.rename(columns={'Positional Score': 'Score', 'DENOVO': 'Sequence'})
    filter_df['Score'] = filter_df['Score'].apply(lambda x: re.sub('[\\[\\]\\s]', '', str(x)))
    filter_df["Score"] = filter_df["Score"].apply(lambda x: list(map(float, [y for y in x.split(',') if y != ''])))
    filter_df = filter_df.dropna(axis=0, how='any')
    filter_df['Sequence'] = filter_df['Sequence'].apply(lambda x: re.sub('[^A-Z]', '', str(x)))

    filter_df["Filtered_Subsequences"] = filter_df.apply(
        lambda row: extract_subsequences_sub(row["Sequence"], row["Score"],
                                            min_score=min_score, min_length=min_length), axis=1)
    filter_df = filter_df[filter_df["Filtered_Subsequences"].str.len() > 0]
    return filter_df

def process_instanovo(path, min_score=0.75, min_length=8):
    files = [x for x in os.listdir(path) if x.endswith('.result.csv')]
    if not files:
        raise FileNotFoundError("No instanovo .result.csv found")
    dfs = []
    for f in files:
        dfs.append(pd.read_csv(os.path.join(path, f)))
    df = pd.concat(dfs, ignore_index=True)

    df = df.rename(columns={'preds': 'DENOVO', 'token_log_probs': 'Positional Score_raw'})
    filter_df = df[['DENOVO', 'Positional Score_raw']].copy()
    filter_df['Positional Score'] = filter_df['Positional Score_raw'].apply(
        lambda x: [float(np.exp(float(val))) for val in str(x)[1:-1].split(',') if val.strip() != ''])
    filter_df = filter_df[['DENOVO', 'Positional Score']]
    filter_df = filter_df.rename(columns={'Positional Score': 'Score', 'DENOVO': 'Sequence'})
    filter_df['Sequence'] = filter_df['Sequence'].apply(lambda x: re.sub('[^A-Z]', '', str(x)))

    filter_df["Filtered_Subsequences"] = filter_df.apply(
        lambda row: extract_subsequences_sub(row["Sequence"], row["Score"],
                                            min_score=min_score, min_length=min_length), axis=1)
    filter_df = filter_df[filter_df["Filtered_Subsequences"].str.len() > 0]
    return filter_df

def process_casanovo(path, min_score=0.75, min_length=8):
    df = pd.DataFrame()
    file_list = [x for x in os.listdir(path) if x.endswith('.mztab')]
    for file in file_list:
        with open(os.path.join(path, file), 'r') as f:
            lines = f.readlines()
        filtered_lines = [line for line in lines if not line.startswith('MTD')]
        filtered_data = StringIO(''.join(filtered_lines))
        df_single = pd.read_csv(filtered_data, sep='\\t')
        df = df_single if df.empty else pd.concat([df, df_single])

    filter_df = df[['sequence', 'opt_ms_run[1]_aa_scores']].copy()
    filter_df['sequence'] = filter_df['sequence'].apply(lambda x: re.sub('[^A-Z]', '', str(x)))
    filter_df = filter_df.rename(columns={'opt_ms_run[1]_aa_scores': 'Score', 'sequence': 'Sequence'})
    filter_df["Score"] = filter_df["Score"].apply(lambda x: list(map(float, str(x).split(','))))

    filter_df["Filtered_Subsequences"] = filter_df.apply(
        lambda row: extract_subsequences_sub(row["Sequence"], row["Score"],
                                            min_score=min_score, min_length=min_length), axis=1)
    filter_df = filter_df[filter_df["Filtered_Subsequences"].str.len() > 0]
    return filter_df

peptide_soft_dic = {}

# PepNet
pepnet_path = 'Pepnet'
if os.path.isdir(pepnet_path) and os.listdir(pepnet_path):
    pepnet_files = [f for f in os.listdir(pepnet_path) if f.endswith('.result')]
    if pepnet_files:
        try:
            df = process_pepnet(pepnet_path, min_score=${params.aa_score_threshold},
                                min_length=${params.min_peptide_length})
            peptides = []
            for _, row in df.iterrows():
                peptides.extend(row['Filtered_Subsequences'])
            for pep in set(peptides):
                peptide_soft_dic[pep] = peptide_soft_dic.get(pep, []) + ['Pepnet']
            print(f"INFO: Pepnet contributed {len(set(peptides))} unique peptides", file=sys.stderr)
        except Exception as e:
            print(f"WARNING: Pepnet processing failed: {e}", file=sys.stderr)

# InstaNovo
instanovo_path = 'Instanovo'
if os.path.isdir(instanovo_path) and os.listdir(instanovo_path):
    try:
        df = process_instanovo(instanovo_path, min_score=${params.aa_score_threshold},
                               min_length=${params.min_peptide_length})
        peptides = []
        for _, row in df.iterrows():
            peptides.extend(row['Filtered_Subsequences'])
        for pep in set(peptides):
            peptide_soft_dic[pep] = peptide_soft_dic.get(pep, []) + ['Instanovo']
        print(f"INFO: InstaNovo contributed {len(set(peptides))} unique peptides", file=sys.stderr)
    except Exception as e:
        print(f"WARNING: InstaNovo processing failed: {e}", file=sys.stderr)

# Casanovo
casanovo_path = 'CasaNovo'
if os.path.isdir(casanovo_path) and os.listdir(casanovo_path):
    casanovo_files = [f for f in os.listdir(casanovo_path) if f.endswith('.mztab')]
    if casanovo_files:
        try:
            df = process_casanovo(casanovo_path, min_score=${params.aa_score_threshold},
                                  min_length=${params.min_peptide_length})
            peptides = []
            for _, row in df.iterrows():
                peptides.extend(row['Filtered_Subsequences'])
            for pep in set(peptides):
                peptide_soft_dic[pep] = peptide_soft_dic.get(pep, []) + ['CasaNovo']
            print(f"INFO: CasaNovo contributed {len(set(peptides))} unique peptides", file=sys.stderr)
        except Exception as e:
            print(f"WARNING: CasaNovo processing failed: {e}", file=sys.stderr)

if len(peptide_soft_dic) == 0:
    print("ERROR: No peptides found from any de novo tool!", file=sys.stderr)
    sys.exit(1)

with open('merged_raw_peptide.fasta', 'w') as f:
    for peptide in peptide_soft_dic.keys():
        f.write(f'>{peptide}\\n{peptide}\\n')

print(f"SUCCESS: Total unique peptides: {len(peptide_soft_dic)}", file=sys.stderr)
PYTHON_EOF
  """
}