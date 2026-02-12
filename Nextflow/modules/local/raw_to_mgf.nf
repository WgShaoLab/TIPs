nextflow.enable.dsl=2

process RAW_TO_MGF {
  tag "${sample_id}:${raw.baseName}"

  input:
  tuple val(sample_id), path(raw)

  output:
  tuple val(sample_id), path("${raw.baseName}.mgf"), emit: mgf

  script:
  """
  set -euo pipefail

  RAW_PATH="\$(realpath ${raw})"
  OUTDIR="\$PWD"

  cd /opt/thermorawfileparser
  mono ThermoRawFileParser.exe \
    -i "\$RAW_PATH" \
    -o "\$OUTDIR" \
    --format=0

  cd "\$OUTDIR"

  if [ ! -f "${raw.baseName}.mgf" ]; then
    f=\$(ls -1 *.mgf 2>/dev/null | head -n 1 || true)
    test -n "\$f" || { echo "[ERROR] No mgf produced"; ls -la; exit 2; }
    mv "\$f" "${raw.baseName}.mgf"
  fi
  """
}