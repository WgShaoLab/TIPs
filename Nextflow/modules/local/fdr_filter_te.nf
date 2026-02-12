nextflow.enable.dsl=2

process FDR_FILTER_TE {
  tag "${sample_id}_${engine}_${type}"

  input:
  tuple val(sample_id), val(engine), val(type), path(merged_pepxml)

  output:
  tuple val(sample_id), val(engine), val(type), path("${engine}_merged_${type}_fdr.csv"), emit: fdr_filtered

  script:
  """
  set -euo pipefail

  /opt/tpp/bin/pepxml2csv.py ${merged_pepxml}

  python << 'PYTHON_EOF'
import sys
sys.path.append('${params.code_base}')
from TE_immunopeptide.TE_searchPipline.SoftPipeline_class_iProphet import filter_by_fdr

csv_file = '${merged_pepxml}'.replace('.xml', '.csv')
filter_by_fdr(csv_file, ${params.fdr_threshold})
PYTHON_EOF
  """
}