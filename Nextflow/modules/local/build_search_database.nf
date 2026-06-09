nextflow.enable.dsl=2

process BUILD_SEARCH_DATABASE {
  tag "${sample_id}_${engine}"

  input:
  tuple val(sample_id), val(engine), path(te_fasta)
  path uniprot_fasta
  path contaminants_fasta

  output:
  tuple val(sample_id), val(engine), path("${engine}_DBsearch.fasta"), emit: database

  script:
  """
  python << 'PYTHON_EOF'
# Match CLI search_te._build_combined_db: pure Python, preserve full headers
from pathlib import Path
from typing import Iterable, Optional

def iter_fasta_records(fa_path: str) -> Iterable[tuple]:
    header: Optional[str] = None
    seq_chunks = []
    with open(fa_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            yield header, "".join(seq_chunks)

def write_fasta(records, out_path):
    with open(out_path, "w", encoding="utf-8") as out:
        for h, s in records:
            out.write(f">{h}\\n{s}\\n")

engine = "${engine}"

# Step 1: concatenate human + TE + contaminants (same order as CLI)
write_fasta(
    (r for fa in ["${uniprot_fasta}", "${te_fasta}", "${contaminants_fasta}"]
       for r in iter_fasta_records(fa)),
    f"{engine}.fasta"
)

# Step 2: append rev_ decoys with reversed sequences
def target_decoy_iter():
    for h, s in iter_fasta_records(f"{engine}.fasta"):
        yield h, s
    for h, s in iter_fasta_records(f"{engine}.fasta"):
        yield f"rev_{h}", s[::-1]

write_fasta(target_decoy_iter(), f"{engine}_DBsearch.fasta")
PYTHON_EOF
  """
}