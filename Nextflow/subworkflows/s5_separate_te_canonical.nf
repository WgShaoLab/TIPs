nextflow.enable.dsl=2

include { SEPARATE_TE_CANONICAL } from '../modules/local/separate_te_canonical'

workflow S5_SEPARATE_TE_CANONICAL {

  take:
  all_pepxml

  main:
  separated = SEPARATE_TE_CANONICAL(all_pepxml)

  emit:
  te_pepxml = separated.te_pepxml
  canonical_pepxml = separated.canonical_pepxml
}