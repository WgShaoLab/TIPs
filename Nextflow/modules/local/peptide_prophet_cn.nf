nextflow.enable.dsl=2

process PEPTIDE_PROPHET_CN {
  tag "${sample_id}_${engine}_${type}"

  input:
  tuple val(sample_id), val(engine), val(type), path(pepxml_files), path(database)

  output:
  tuple val(sample_id), val(engine), val(type), path("${engine}_merged_${type}.pep.xml"), emit: merged_pepxml

  script:
  """
  /opt/tpp/bin/xinteract \\
    -OAPNd -PPM -eN -p0 \\
    -THREADS=${params.threads} \\
    -drev_ \\
    -D${database} \\
    -N${engine}_merged_${type}.pep.xml \\
    ${pepxml_files}
  """
}