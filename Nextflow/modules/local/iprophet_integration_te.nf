nextflow.enable.dsl=2

process IPROPHET_INTEGRATION_TE {
  tag "${sample_id}_${type}"

  input:
  tuple val(sample_id), val(type), path(comet_pepxml), path(msfragger_pepxml), path(msgf_pepxml), path(database)

  output:
  tuple val(sample_id), val(type), path("interact.iproph.pep.xml"), emit: iprophet_pepxml
  tuple val(sample_id), path("peptide.tsv"), emit: peptide_tsv

  script:
  """
  set -euo pipefail
  
  philosopher workspace --clean
  philosopher workspace --init
  philosopher database --annotate ${database}

  /opt/tpp/bin/RefreshParser ${comet_pepxml} ${database}
  /opt/tpp/bin/RefreshParser ${msfragger_pepxml} ${database}
  /opt/tpp/bin/RefreshParser ${msgf_pepxml} ${database}

  philosopher iprophet --threads ${params.threads} \\
    ${comet_pepxml} ${msfragger_pepxml} ${msgf_pepxml}

  philosopher proteinprophet --iproph.pep.xml interact.iproph.pep.xml

  philosopher filter --pepxml interact.iproph.pep.xml --pep ${params.fdr_threshold} --psm ${params.fdr_threshold}
  philosopher report
  """
}