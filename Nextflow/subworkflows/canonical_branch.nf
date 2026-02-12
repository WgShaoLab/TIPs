nextflow.enable.dsl=2

include { PEPTIDE_PROPHET_CN }       from '../modules/local/peptide_prophet_cn'
include { FDR_FILTER_CN }            from '../modules/local/fdr_filter_cn'
include { IPROPHET_INTEGRATION_CN }  from '../modules/local/iprophet_integration_cn'

workflow CANONICAL_BRANCH {

  take:
  canonical_pepxml
  search_databases

  main:
  canonical_prophet_input = canonical_pepxml
    .map { sid, engine, files -> tuple(sid, engine, 'Canonical', files) }
    .combine(
      search_databases.map { sid, engine, db -> tuple(sid, engine, db) },
      by: [0, 1]
    )
    .map { sid, engine, type, files, db -> tuple(sid, engine, type, files, db) }

  canonical_prophet = PEPTIDE_PROPHET_CN(canonical_prophet_input)

  canonical_fdr = FDR_FILTER_CN(canonical_prophet.merged_pepxml)

  canonical_iprophet_input = canonical_prophet.merged_pepxml
    .map { sid, engine, type, pepxml -> tuple(sid, type, engine, pepxml) }
    .groupTuple(by: [0, 1])
    .map { sid, type, engines2, pepxmls ->
      def engineMap = [:]
      engines2.eachWithIndex { engine, idx -> engineMap[engine] = pepxmls[idx] }

      if (!engineMap.containsKey('Comet') ||
          !engineMap.containsKey('MSfragger') ||
          !engineMap.containsKey('MS_GF')) {
        error "Missing engine output for ${sid}_${type}: found ${engines2}"
      }

      tuple(sid, type, engineMap['Comet'], engineMap['MSfragger'], engineMap['MS_GF'])
    }
    .combine(
      search_databases
        .filter { sid, engine, db -> engine == 'MS_GF' }
        .map { sid, engine, db -> tuple(sid, db) },
      by: 0
    )
    .map { sid, type, comet_xml, msfragger_xml, msgf_xml, db ->
      tuple(sid, type, comet_xml, msfragger_xml, msgf_xml, db)
    }

  canonical_iprophet = IPROPHET_INTEGRATION_CN(canonical_iprophet_input)

  emit:
  iprophet_pepxml = canonical_iprophet.iprophet_pepxml
  peptide_tsv = canonical_iprophet.peptide_tsv
  canonical_fdr = canonical_fdr.fdr_filtered
}