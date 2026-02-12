nextflow.enable.dsl=2

include { BLASTP_DENOVO_PEPTIDES } from '../modules/local/blastp_denovo_peptides'
include { BUILD_DENOVO_DATABASE }  from '../modules/local/build_denovo_database'

workflow S2_BLASTP_BUILD_TE_DB {

  take:
  peptides
  te_fasta
  te_npy
  te_class_dic

  main:
  blastp_results = BLASTP_DENOVO_PEPTIDES(peptides, te_fasta)

  denovo_database = BUILD_DENOVO_DATABASE(
    blastp_results.blastp_dir.join(peptides),
    te_npy,
    te_class_dic
  )

  emit:
  te_fasta = denovo_database.te_fasta
  tagged_fasta = denovo_database.tagged_fasta
}