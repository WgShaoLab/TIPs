nextflow.enable.dsl=2

process MSFRAGGER_SEARCH {
  tag "${sample_id}"
  maxForks 1

  input:
  tuple val(sample_id), path(mzml_files), path(database)
  path msfragger_params_template

  output:
  tuple val(sample_id), val('MSfragger'), path("*.pepXML"), emit: pepxml

  script:
  """
  set -euo pipefail
  sed "s|database_name =.*|database_name = ${database}|" \\
    ${msfragger_params_template} > msfragger.params

  java -Xmx120g -jar /opt/msfragger/MSFragger-4.1.jar \\
    msfragger.params ${mzml_files}
  """
}