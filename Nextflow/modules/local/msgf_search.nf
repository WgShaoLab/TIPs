nextflow.enable.dsl=2

process MSGF_SEARCH {
  tag "${sample_id}_${mzml.baseName}"
  maxForks 1

  input:
  tuple val(sample_id), path(mzml), path(database)
  path msgf_params

  output:
  tuple val(sample_id), path("*.mzid"), emit: mzid

  script:
  """
  set -euo pipefail
  java -Xmx120g -jar /opt/msgf/MSGFPlus.jar \\
    -s ${mzml} \\
    -conf "${msgf_params}" \\
    -d "${database}" \\
    -o "${mzml.baseName}.mzid"
  """
}