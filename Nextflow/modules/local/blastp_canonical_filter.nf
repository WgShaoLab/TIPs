nextflow.enable.dsl=2

process BLASTP_CANONICAL_FILTER {
    tag "${sid}_blastp_filter"

    input:
    tuple val(sid), path(peptide_tsv)
    path blastp_db

    output:
    tuple val(sid), path("peptide_blastp.filter.txt"), emit: filtered_peptides

    script:
    """
    #!/bin/bash
    set -euo pipefail

    # Find the actual BLAST DB prefix by globbing for index files
    export DB_PREFIX=\$(ls "${blastp_db}"/*.pin 2>/dev/null | head -1 | sed 's/\\.[0-9][0-9]\\.pin\$//; s/\\.pin\$//')
    if [ -z "\${DB_PREFIX}" ]; then
      echo "WARNING: No BLAST DB (.pin) found in ${blastp_db}, trying alternative patterns..."
      export DB_PREFIX=\$(ls "${blastp_db}"/*.phr 2>/dev/null | head -1 | sed 's/\\.[0-9][0-9]\\.phr\$//; s/\\.phr\$//')
    fi
    if [ -z "\${DB_PREFIX}" ]; then
      echo "ERROR: Cannot locate BLAST database in ${blastp_db}"
      ls -la "${blastp_db}/" || true
      exit 1
    fi
    echo "Using BLAST database: \${DB_PREFIX}"

    python3 << 'PYEOF'
import pandas as pd
import subprocess
import os

df = pd.read_csv("${peptide_tsv}", sep='\\t')
peptide_col = 'peptide' if 'peptide' in df.columns else ('Peptide' if 'Peptide' in df.columns else None)

if not peptide_col or df.empty:
    df.to_csv("peptide_blastp.filter.txt", sep='\\t', index=False)
    print("No peptide column found or empty dataframe, skipping BLASTP")
    exit(0)

peptides = list(set(df[peptide_col].dropna().tolist()))
if not peptides:
    df.to_csv("peptide_blastp.filter.txt", sep='\\t', index=False)
    print("No peptides found, skipping BLASTP")
    exit(0)

print(f"Total unique peptides to check: {len(peptides)}")

with open("query.fasta", 'w') as f:
    for seq in peptides:
        f.write(f">{seq}\\n{seq}\\n")

db_path = os.environ.get("DB_PREFIX", "")
if not db_path:
    print("ERROR: DB_PREFIX not set", file=__import__('sys').stderr)
    exit(1)
print(f"Using BLAST database: {db_path}")

cmd = [
    'blastp', '-task', 'blastp-short',
    '-query', 'query.fasta',
    '-db', db_path,
    '-out', 'blastp_result.txt',
    '-outfmt', '6',
    '-evalue', '20000',
    '-num_threads', '${task.cpus ?: 8}'
]

print(f"Running: {' '.join(cmd)}")
result = subprocess.run(cmd, capture_output=True, text=True)
if result.returncode != 0:
    print(f"BLASTP stderr: {result.stderr}")

remove_set = set()
if os.path.exists("blastp_result.txt") and os.path.getsize("blastp_result.txt") > 0:
    cols = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']
    bp = pd.read_csv("blastp_result.txt", sep='\\t', names=cols)
    bp['qlen'] = bp['qseqid'].apply(len)

    bp_match = bp[(bp['qstart'] == 1) & (bp['qend'] == bp['qlen']) & (bp['pident'] == 100)]
    remove_set = set(bp_match['qseqid'].tolist())
    print(f"Found {len(remove_set)} peptides matching canonical proteins")

df_out = df[~df[peptide_col].isin(remove_set)]
df_out.to_csv("peptide_blastp.filter.txt", sep='\\t', index=False)
print(f"Removed {len(remove_set)} canonical peptides, kept {len(df_out)}/{len(df)} rows")
PYEOF
    """
}