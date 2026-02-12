nextflow.enable.dsl=2

process DENOVO_ALL_BATCH {
  tag "${sample_id}"

  input:
  tuple val(sample_id), path(mgf_files)

  output:
  tuple val(sample_id), path("casanovo/*.mztab"), emit: casanovo_results
  tuple val(sample_id), path("pepnet/*.result"), emit: pepnet_results
  tuple val(sample_id), path("instanovo/*.result.csv"), emit: instanovo_results

  script:
  def mgf_args = mgf_files.collect { "\"${it}\"" }.join(' ')

  """
  set -euo pipefail

  test -f "${params.denovo_runner}" || { echo "[ERROR] denovo_runner not found: ${params.denovo_runner}" >&2; exit 2; }
  python3 -V

  MGF_DIR=\$(python3 - ${mgf_args} <<'PY'
from pathlib import Path
import sys

mgfs = [Path(p).resolve() for p in sys.argv[1:]]
if not mgfs:
    raise SystemExit("[ERROR] empty mgf list")

parents = {p.parent for p in mgfs}
if len(parents) != 1:
    raise SystemExit(f"[ERROR] mgf files are not in the same directory: {parents}")

print(str(next(iter(parents))))
PY
  )

  echo "[INFO] sample_id=${sample_id}"
  echo "[INFO] MGF_DIR=\$MGF_DIR"

  mkdir -p denovo_out

  python3 "${params.denovo_runner}" \\
    --tool all \\
    --input "\$MGF_DIR" \\
    --outdir denovo_out \\
    --sample "${sample_id}" \\
    ##### --param-file "${params.casanovo_config}"

  mkdir -p casanovo pepnet instanovo
  cp -a denovo_out/casanovo/*.mztab casanovo/ 2>/dev/null || true
  cp -a denovo_out/pepnet/*.result pepnet/ 2>/dev/null || true
  cp -a denovo_out/instanovo/*.result.csv instanovo/ 2>/dev/null || true

  echo "[INFO] denovo outputs:"
  ls -lh casanovo || true
  ls -lh pepnet || true
  ls -lh instanovo || true
  """
}