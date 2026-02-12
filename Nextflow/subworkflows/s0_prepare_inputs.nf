nextflow.enable.dsl=2

include { RAW_TO_MZML } from '../modules/local/raw_to_mzml'
include { RAW_TO_MGF }  from '../modules/local/raw_to_mgf'

workflow S0_PREPARE_INPUTS {

  take:
  sample_dir

  main:

  sample_id_ch = Channel.value(sample_dir.getName())

  raw_ch = Channel
    .fromPath("${sample_dir}/raw/*.raw", checkIfExists: true)
    .combine(sample_id_ch)
    .map { raw, sid -> tuple(sid, raw) }

  mzml_each = RAW_TO_MZML(raw_ch)
  mgf_each  = RAW_TO_MGF(raw_ch)

  // Build sorted lists (single tuple per sample_id)
  mzml_list = mzml_each.mzml
    .map { sid, mzml -> tuple(sid, mzml) }
    .groupTuple()
    .map { sid, files -> tuple(sid, files.sort()) }

  mgf_list = mgf_each.mgf
    .map { sid, mgf -> tuple(sid, mgf) }
    .groupTuple()
    .map { sid, files -> tuple(sid, files.sort()) }

  emit:
  sample_id  = sample_id_ch
  mzml_files = mzml_list
  mgf_files  = mgf_list
}