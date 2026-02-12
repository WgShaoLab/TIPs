nextflow.enable.dsl=2

process BUILD_SEARCH_DATABASE {
  tag "${sample_id}_${engine}"

  input:
  tuple val(sample_id), val(engine), path(te_fasta_tagged)
  path uniprot_fasta
  path contaminants_fasta

  output:
  tuple val(sample_id), val(engine), path("${engine}_DBsearch.fasta"), emit: database

  script:
  """
  set -euo pipefail

  cat ${uniprot_fasta} ${te_fasta_tagged} ${contaminants_fasta} > ${engine}.fasta

  python << 'PYTHON_EOF'
from Bio import SeqIO

with open('${engine}_decoy.fasta', 'w') as out:
    for record in SeqIO.parse('${engine}.fasta', 'fasta'):
        decoy_seq = str(record.seq)[::-1]
        out.write(f'>rev_{record.id}\\n{decoy_seq}\\n')
PYTHON_EOF

  cat ${engine}.fasta ${engine}_decoy.fasta > ${engine}_DBsearch.fasta
  """
}