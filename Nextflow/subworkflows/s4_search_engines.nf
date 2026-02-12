nextflow.enable.dsl=2

include { COMET_SEARCH }      from '../modules/local/comet_search'
include { MSFRAGGER_SEARCH }  from '../modules/local/msfragger_search'
include { MSGF_SEARCH }       from '../modules/local/msgf_search'
include { MSGF_CONVERT }      from '../modules/local/msgf_convert'

workflow S4_SEARCH_ENGINES {

  take:
  mzml_files
  search_databases
  comet_params
  msfragger_params
  msgf_params

  main:
  // Comet
  comet_db = search_databases
    .filter { sid, engine, db -> engine == 'Comet' }
    .map { sid, engine, db -> tuple(sid, db) }

  comet_input = mzml_files
    .combine(comet_db, by: 0)  // by sample_id
    .map { sid, files, db -> tuple(sid, files, db) }

  comet_pepxml = COMET_SEARCH(comet_input, comet_params)

  // MSFragger
  msfragger_db = search_databases
    .filter { sid, engine, db -> engine == 'MSfragger' }
    .map { sid, engine, db -> tuple(sid, db) }

  msfragger_input = mzml_files
    .combine(msfragger_db, by: 0)
    .map { sid, files, db -> tuple(sid, files, db) }

  msfragger_pepxml = MSFRAGGER_SEARCH(msfragger_input, msfragger_params)

  // MSGF (per mzML file)
  msgf_db = search_databases
    .filter { sid, engine, db -> engine == 'MS_GF' }
    .map { sid, engine, db -> tuple(sid, db) }

  msgf_input = mzml_files
    .combine(msgf_db, by: 0)
    .flatMap { sid, files, db ->
      files.collect { mzml_file -> tuple(sid, mzml_file, db) }
    }

  msgf_mzid = MSGF_SEARCH(msgf_input, msgf_params)

  msgf_pepxml = MSGF_CONVERT(msgf_mzid.mzid.groupTuple())

  all_pepxml = comet_pepxml.pepxml
    .mix(msfragger_pepxml.pepxml, msgf_pepxml.pepxml)

  emit:
  all_pepxml = all_pepxml
}