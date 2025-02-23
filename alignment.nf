#!/usr/bin/env nextflow
nextflow.enable.dsl=2

import groovy.json.JsonSlurper

include { INDEX_GENOME             } from './subworkflows/index_genome/main'
include { FASTQ_FASTP_FASTQC       } from './subworkflows/fastq_fastp_fastqc/main'
include { ALIGN_MARKDUP_BQSR_STATS } from './subworkflows/align_markdup_bqsr_stats/main'

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

    // If params.input_json exists, use that for multi sample processing.
    // Otherwise, default to single sample analysis.
    multi_params = params.input_json ? jsonSlurper.parse(new File(params.input_json)).collect{params + it} : [params]
    output_dir = params.output_dir ? params.output_dir : "."
    
    // Process samples through the pipeline
    samples = Channel.from(multi_params.collect{ it -> tuple([
                id: it.specimen_id, pid: it.patient_id, single_end:false, tissue: it.tissue, purity: it.purity ],
                [ file(it.read1, checkIfExists: true), file(it.read2, checkIfExists: true) ]) })
    ch_versions = Channel.empty()

    // Index the reference genome for alignment with BWA mem
    INDEX_GENOME ( [[ id:'genome'], file(params.reference_file, checkIfExists: true)])

    // Trim sequence reads with paired-end data
    FASTQ_FASTP_FASTQC ( samples,
                        file(params.adapter_fasta, checkIfExists:true),
                        params.save_trimmed_fail,
                        params.save_merged,
                        params.skip_fastp,
                        params.skip_fastqc )
    ch_versions = ch_versions.mix( FASTQ_FASTP_FASTQC.out.versions )

    // Align reads to a reference genome and generate alignment statistics
    ALIGN_MARKDUP_BQSR_STATS ( FASTQ_FASTP_FASTQC.out.reads,
                INDEX_GENOME.out.index,
                file(params.reference_file, checkIfExists: true),
                params.val_sort_bam,
                file(params.exome_plus_tumor_panel_bed, checkIfExists: true),
                INDEX_GENOME.out.fai,
                INDEX_GENOME.out.dict,
                file(params.known_snp_vcf, checkIfExists: true),
                file(params.known_snp_vcf_tbi, checkIfExists: true),
                file(params.known_indel_vcf, checkIfExists: true),
                file(params.known_indel_vcf_tbi, checkIfExists: true),
                params.split_reads)
    ch_versions = ch_versions.mix( ALIGN_MARKDUP_BQSR_STATS.out.versions )
}