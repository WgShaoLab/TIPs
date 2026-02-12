nextflow.enable.dsl=2

process COMET_SEARCH {
  tag "${sample_id}"
  maxForks 1

  input:
  tuple val(sample_id), path(mzml_files), path(database)
  path comet_params

  output:
  tuple val(sample_id), val('Comet'), path("*.pep.xml"), emit: pepxml

  script:
  """
  set -euo pipefail
  /opt/tpp/bin/comet -P${comet_params} -D${database} ${mzml_files}
  """
}