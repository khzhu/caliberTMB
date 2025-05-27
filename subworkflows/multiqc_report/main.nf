#!/usr/bin/env nextflow

include { MULTIQC } from '../../modules/multiqc/main'

workflow MULTIQC_REPORT {

    take:
    ch_bam_files
    ch_input_qc_files

    main:
    ch_versions         = Channel.empty()

    MULTIQC ( ch_input_qc_files )
    ch_versions = ch_versions.mix(MULTIQC.out.versions)

    emit:
    report   = MULTIQC.out.report  // channel: [ val(meta), path(report) ]
    data     = MULTIQC.out.data    // channel: [ val(meta), path(data) ]
    versions = ch_versions         // channel: [ path(versions.yml) ]
}