nextflow.enable.dsl=2

include { PEPTIDE_PROPHET_TE }       from '../modules/local/peptide_prophet_te'
include { FDR_FILTER_TE }            from '../modules/local/fdr_filter_te'
include { IPROPHET_INTEGRATION_TE }  from '../modules/local/iprophet_integration_te'
include { BLASTP_CANONICAL_FILTER }  from '../modules/local/blastp_canonical_filter'
include { MHC_BINDING_PREDICTION }   from '../modules/local/mhc_binding_prediction'

workflow TE_BRANCH {

  take:
  te_pepxml
  search_databases
  blastp_db

  main:
  // PeptideProphet TE
  te_prophet_input = te_pepxml
    .map { sid, engine, files -> tuple(sid, engine, 'TE', files) }
    .combine(
      search_databases.map { sid, engine, db -> tuple(sid, engine, db) },
      by: [0, 1]
    )
    .map { sid, engine, type, files, db -> tuple(sid, engine, type, files, db) }

  te_prophet = PEPTIDE_PROPHET_TE(te_prophet_input)

  // FDR filter (kept as QC output; not used downstream)
  te_fdr = FDR_FILTER_TE(te_prophet.merged_pepxml)

  // iProphet TE
  te_iprophet_input = te_prophet.merged_pepxml
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

  te_iprophet = IPROPHET_INTEGRATION_TE(te_iprophet_input)

  // BLASTP filter + MHC prediction
  filtered = BLASTP_CANONICAL_FILTER(
    te_iprophet.peptide_tsv.map { sid, tsv -> tuple(sid, tsv) },
    blastp_db
  )

  final_results = MHC_BINDING_PREDICTION(filtered.filtered_peptides)

  emit:
  final_results = final_results.final_results
  iprophet_pepxml = te_iprophet.iprophet_pepxml
  peptide_tsv = te_iprophet.peptide_tsv
  te_fdr = te_fdr.fdr_filtered
}