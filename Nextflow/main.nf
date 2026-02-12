#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { TIPS } from './workflows/tips'

if (params.help) {
  log.info """
  TIPs modular pipeline (nf-core style layout)

  Usage:
    nextflow run main.nf --sample_dir <DIR> --outdir <DIR>

  Required:
    --sample_dir   Directory containing raw/*.raw

  Outputs:
    --outdir       Output directory (default: ${params.outdir})
  """.stripIndent()
  exit 0
}

if (!params.sample_dir) {
  error "Please provide --sample_dir (a directory containing raw/*.raw)"
}

workflow {
  TIPS( file(params.sample_dir) )
}