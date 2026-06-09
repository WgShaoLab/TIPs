nextflow.enable.dsl=2

include { SEPARATE_TE_CANONICAL } from '../modules/local/separate_te_canonical'

workflow S5_SEPARATE_TE_CANONICAL {

  take:
  all_pepxml     // tuple(sample_id, engine, pepxml_files)
  te_fasta_ch    // tuple(sample_id, te_fasta_path)

  main:
  // Join TE FASTA path with pepXML by sample_id so each sample gets its own TE FASTA
  separated_input = all_pepxml
    .combine(te_fasta_ch, by: 0)
    .map { sid, engine, files, te_fa -> tuple(sid, engine, files, te_fa) }

  separated = SEPARATE_TE_CANONICAL(separated_input)

  emit:
  te_pepxml = separated.te_pepxml
  canonical_pepxml = separated.canonical_pepxml
}