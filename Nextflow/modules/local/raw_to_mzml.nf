nextflow.enable.dsl=2

process RAW_TO_MZML {
  tag "${sample_id}:${raw.baseName}"

  input:
  tuple val(sample_id), path(raw)

  output:
  tuple val(sample_id), path("${raw.baseName}.mzML"), emit: mzml

  script:
  """
  set -euo pipefail

  RAW_PATH="\$(realpath ${raw})"
  OUTDIR="\$PWD"

  cd /opt/thermorawfileparser
  mono ThermoRawFileParser.exe \
    -i "\$RAW_PATH" \
    -o "\$OUTDIR" \
    --format=1

  cd "\$OUTDIR"

  if [ ! -f "${raw.baseName}.mzML" ]; then
    f=\$(ls -1 *.mzML *.mzml 2>/dev/null | head -n 1 || true)
    test -n "\$f" || { echo "[ERROR] No mzML produced"; ls -la; exit 2; }
    mv "\$f" "${raw.baseName}.mzML"
  fi
  """
}