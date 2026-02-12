nextflow.enable.dsl=2

include { BUILD_SEARCH_DATABASE } from '../modules/local/build_search_database'

workflow S3_BUILD_SEARCH_DB {

  take:
  te_fasta_tagged
  uniprot_fasta
  contaminants_fasta

  main:
  engines = Channel.of('Comet', 'MSfragger', 'MS_GF')

  // te_fasta_tagged is: tuple(sid, te_fasta)
  search_databases = BUILD_SEARCH_DATABASE(
    engines
      .combine(te_fasta_tagged)
      .map { engine, sid, te_fasta -> tuple(sid, engine, te_fasta) },
    uniprot_fasta,
    contaminants_fasta
  )

  emit:
  search_databases = search_databases.database
}