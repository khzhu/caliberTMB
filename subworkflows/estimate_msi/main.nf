#!/usr/bin/env nextflow

include { MSISENSOR2_MSI } from '../../modules/msisensor2/main'

workflow ESTIMATE_MSI {

    take:
    tumor_bam_ch    // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    models          // channel: [mandatory] [ path(models) ]

    main:
    ch_versions         = Channel.empty()

    MSISENSOR2_MSI ( tumor_bam_ch, models )
    ch_versions = ch_versions.mix(MSISENSOR2_MSI.out.versions)

    emit:
    msi          = MSISENSOR2_MSI.out.msi          // channel: [ val(meta), path(msi) ]
    distribution = MSISENSOR2_MSI.out.distribution // channel: [ val(meta), path(distribution) ]
    somatic      = MSISENSOR2_MSI.out.somatic      // channel: [ val(meta), path(somatic) ]
    versions     = ch_versions                     // channel: [ path(versions.yml) ]
}