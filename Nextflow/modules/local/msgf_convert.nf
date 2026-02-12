nextflow.enable.dsl=2

process MSGF_CONVERT {
  tag "${sample_id}"

  input:
  tuple val(sample_id), path(mzid_files)

  output:
  tuple val(sample_id), val('MS_GF'), path("*.pep.xml"), emit: pepxml

  script:
  """
  set -euo pipefail

  for mzid in ${mzid_files}; do
    base=\$(basename \$mzid .mzid)

    /opt/scripts/mzid_conver.sh \$mzid \${base}

    /opt/pwiz/idconvert \${base}_modify.mzid --pepXML -o .

    if [ -f "\${base}.pepXML" ]; then
      mv "\${base}.pepXML" "\${base}.pep.xml"
    else
      echo "ERROR: Expected file \${base}.pepXML not found!" >&2
      ls -la >&2
      exit 1
    fi
  done

  ls -lh *.pep.xml >&2
  """
}