#!/usr/bin/env nextflow
nextflow.enable.dsl=2

import groovy.json.JsonSlurper

include { SNV_MUTECT2  } from './subworkflows/snv_mutect2/main'
include { SNV_STRELKA2 } from './subworkflows/snv_strelka2/main'
include { TMB_CALIBER  } from './subworkflows/calculate_tmb/main'

// Main workflow
workflow {
    log.info """\
    TMB Estimation Pipeline ${params.release}
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    sample sheet  : ${params.input_json}
    output dir    : ${params.output_dir}
    """

    // Read the input sample configuration file
    def jsonSlurper = new JsonSlurper()

    // If params.input_json is available, use it for multi-sample processing.
    // If not, default to single sample analysis.
    multi_params = params.input_json ? jsonSlurper.parse(new File(params.input_json)).collect{params + it} : [params]
    output_dir = params.output_dir ? params.output_dir : "."
    
    // Process samples through the pipeline
    samples = Channel.from(multi_params.collect{ it -> tuple([
                id: it.patient_id, tissue: it.tissue, purity: it.purity ],
                [ file(it.bam_tumor, checkIfExists: true), file(it.bam_normal, checkIfExists: true) ],
                [ file(it.bai_tumor, checkIfExists: true), file(it.bai_normal, checkIfExists: true) ] ) })
    ch_versions = Channel.empty()

    bed_files = Channel.fromPath(params.tumor_panel_bed_files, checkIfExists: true)
    ch_input_files = samples.combine(bed_files)
                    .map { meta, input_bams, input_index_files, intervals ->
                        new_meta = meta.clone()
                        new_meta.sid = intervals.baseName != "no_intervals" ? new_meta.id + "_" + intervals.baseName : new_meta.id
                        intervals = intervals.baseName != "no_intervals" ? intervals : []
                        [new_meta, input_bams, input_index_files, intervals]
                    }
    // Somatic variant calling with Mutect2
    SNV_MUTECT2 (ch_input_files,
                [[ id:'genome'], file(params.reference_file, checkIfExists: true)],
                [[ id:'genome'], file(params.fai_file, checkIfExists: true)],
                [[ id:'genome'], file(params.dict_file, checkIfExists: true)],
                [[id:'gnomad'],file(params.gnomad_exome_vcf, checkIfExists: true),
                    file(params.gnomad_exome_vcf_tbi, checkIfExists: true)],
                [[id:'exac'], file(params.exac_common_vcf, checkIfExists: true),
                    file(params.exac_common_vcf_tbi, checkIfExists: true)],
                [[id: 'cosmic'], file(params.cosmic_vcf, checkIfExists: true),
                    file(params.cosmic_vcf_tbi, checkIfExists: true)],
                file(params.vep_cache, checkIfExists: true))
    ch_versions = ch_versions.mix( SNV_MUTECT2.out.versions )

    // Somatic variant calling with Strelka2
    SNV_STRELKA2 (samples,
                [[ id:'genome'], file(params.reference_file, checkIfExists: true)],
                [[ id:'genome'], file(params.fai_file, checkIfExists: true)],
                [[ id:'genome'], file(params.dict_file, checkIfExists: true)],
                [[ id:'genome'], file(params.tumor_panel_bed_gz, checkIfExists: true),
                                    file(params.tumor_panel_bed_tbi, checkIfExists: true)],
                [[id: 'gnomad'], file(params.gnomad_exome_vcf, checkIfExists: true),
                                    file(params.gnomad_exome_vcf_tbi, checkIfExists: true)],
                [[id: 'cosmic'], file(params.cosmic_vcf, checkIfExists: true),
                                    file(params.cosmic_vcf_tbi, checkIfExists: true)],
                file(params.vep_cache, checkIfExists: true))
    ch_versions = ch_versions.mix( SNV_STRELKA2.out.versions )

    // Gather annotated variants in MAF format files
    mutect2_maf = SNV_MUTECT2.out.maf
        .map { it -> it[1] }
        .collectFile( name: 'mutect2_merged.maf', keepHeader:true, skip:2, storeDir:params.store_dir )
        .map { [ [ id:'tmb'], it ] }

    strelka2_maf = SNV_STRELKA2.out.maf
        .map { it -> it[1] }
        .collectFile( name: 'strelka2_merged.maf', keepHeader:false, skip:2, storeDir:params.store_dir )
        .map { [ [ id:'tmb'], it ] }

    // Estimating tumor mutation burden (TMB)
    TMB_CALIBER ( mutect2_maf, strelka2_maf )
    ch_versions = ch_versions.mix(TMB_CALIBER.out.versions)
}