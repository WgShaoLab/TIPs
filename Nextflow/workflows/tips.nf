nextflow.enable.dsl=2

include { S0_PREPARE_INPUTS }        from '../subworkflows/s0_prepare_inputs'
include { S1_DENOVO }               from '../subworkflows/s1_denovo'
include { S2_BLASTP_BUILD_TE_DB }   from '../subworkflows/s2_blastp_build_te_db'
include { S3_BUILD_SEARCH_DB }      from '../subworkflows/s3_build_search_db'
include { S4_SEARCH_ENGINES }       from '../subworkflows/s4_search_engines'
include { S5_SEPARATE_TE_CANONICAL }from '../subworkflows/s5_separate_te_canonical'
include { TE_BRANCH }               from '../subworkflows/te_branch'
include { CANONICAL_BRANCH }        from '../subworkflows/canonical_branch'

workflow TIPS {

  take:
  sample_dir

  main:
  // Reference channels (all paths are configured in nextflow.config params)
  te_blastdb_ch    = Channel.value(file(params.te_blastdb_dir, type: 'dir', checkIfExists: true))
  te_npy_ch        = Channel.value(file(params.te_npy, checkIfExists: true))
  te_class_dic_ch  = Channel.value(file(params.te_class_dic, checkIfExists: true))
  uniprot_ch       = Channel.value(file(params.uniprot_fasta, checkIfExists: true))
  contaminants_ch  = Channel.value(file(params.contaminants_fasta, checkIfExists: true))
  blastp_db_ch     = Channel.value(file(params.blastp_db))
  comet_params_ch  = Channel.value(file(params.comet_params, checkIfExists: true))
  msfragger_params_ch = Channel.value(file(params.msfragger_params, checkIfExists: true))
  msgf_params_ch   = Channel.value(file(params.msgf_params, checkIfExists: true))

  // S0: raw -> mzML + mgf
  prepared = S0_PREPARE_INPUTS(sample_dir)

  // S1: de novo
  denovo = S1_DENOVO(
    prepared.mgf_files,
    te_npy_ch,
    te_class_dic_ch
  )

  // S2: BLASTP + TE DB build
  te_db = S2_BLASTP_BUILD_TE_DB(
    denovo.peptides,
    te_blastdb_ch,
    te_npy_ch,
    te_class_dic_ch
  )

  // S3: build search DB for engines
  searchdb = S3_BUILD_SEARCH_DB(
    te_db.tagged_fasta,
    uniprot_ch,
    contaminants_ch
  )

  // S4: run search engines (Comet/MSFragger/MSGF)
  searches = S4_SEARCH_ENGINES(
    prepared.mzml_files,
    searchdb.search_databases,
    comet_params_ch,
    msfragger_params_ch,
    msgf_params_ch
  )

  // S5: separate TE and Canonical pepxml
  separated = S5_SEPARATE_TE_CANONICAL(searches.all_pepxml)

  // Branches (from S5)
  te_out = TE_BRANCH(
    separated.te_pepxml,
    searchdb.search_databases,
    blastp_db_ch
  )

  //cn_out = CANONICAL_BRANCH(
  //  separated.canonical_pepxml,
  //  searchdb.search_databases
  //)

  // Show final TE output
  te_out.final_results.view { sid, f ->
    "Sample ${sid}: Final TE peptides saved to ${f}"
  }

  emit:
  final_results = te_out.final_results
  iprophet_te_pepxml = te_out.iprophet_pepxml
  peptide_tsv_te = te_out.peptide_tsv
  //iprophet_cn_pepxml = cn_out.iprophet_pepxml
  //peptide_tsv_cn = cn_out.peptide_tsv
}

workflow.onComplete {
  log.info """
  ===============================================
  TIPs modular pipeline completed!
  ===============================================
  Status:    ${workflow.success ? 'SUCCESS' : 'FAILED'}
  Exit code: ${workflow.exitStatus}
  Duration:  ${workflow.duration}
  Sample:    ${params.sample_dir}
  Outdir:    ${params.outdir}
  ===============================================
  """.stripIndent()
}