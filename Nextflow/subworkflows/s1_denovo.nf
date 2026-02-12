nextflow.enable.dsl=2

include { DENOVO_ALL_BATCH }        from '../modules/local/denovo_all_batch'
include { PROCESS_DENOVO_RESULTS }  from '../modules/local/process_denovo_results'

workflow S1_DENOVO {

  take:
  mgf_files
  te_npy
  te_class_dic

  main:
  denovo_all = DENOVO_ALL_BATCH(mgf_files)

  denovo_combined = denovo_all.casanovo_results
    .join(denovo_all.pepnet_results)
    .join(denovo_all.instanovo_results)

  processed = PROCESS_DENOVO_RESULTS(
    denovo_combined,
    te_npy,
    te_class_dic
  )

  emit:
  peptides = processed.peptides
}