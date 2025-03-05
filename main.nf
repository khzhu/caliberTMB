#!/usr/bin/env nextflow
nextflow.enable.dsl=2

import groovy.json.JsonSlurper

include { INDEX_GENOME             } from './subworkflows/index_genome/main'
include { FASTQ_FASTP_FASTQC       } from './subworkflows/fastq_fastp_fastqc/main'
include { ALIGN_MARKDUP_BQSR_STATS } from './subworkflows/align_markdup_bqsr_stats/main'
include { ESTIMATE_MSI             } from './subworkflows/estimate_msi/main'
include { SNV_MUTECT2              } from './subworkflows/snv_mutect2/main'
include { SNV_STRELKA2             } from './subworkflows/snv_strelka2/main'
include { TMB_CALIBER              } from './subworkflows/calculate_tmb/main'
include { MERGE_MSI                } from './modules/collectmsi/main'

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

    // Provide a JSON file as an input parameter for processing multiple samples
    // If params.input_json is available, use it for multi-sample processing.
    // If not, default to single sample analysis.
    multi_params = params.input_json ? jsonSlurper.parse(new File(params.input_json)).collect{params + it} : [params]
    output_dir = params.output_dir ? params.output_dir : "."
    
    // Process samples through the pipeline
    samples = Channel.from(multi_params.collect{ it -> tuple([
                id: it.specimen_id, pid: it.patient_id, single_end:false, tissue: it.tissue, purity: it.purity ],
                [ file(it.read1, checkIfExists: true), file(it.read2, checkIfExists: true) ]) })
    ch_versions = Channel.empty()

    // Index the reference genome for alignment with BWA mem
    INDEX_GENOME ( [[ id:'genome'], file(params.reference_file, checkIfExists: true)])

    // Trim paired-end sequence reads and generate a QC report
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
                file(params.tumor_panel_bed, checkIfExists: true),
                INDEX_GENOME.out.fai,
                INDEX_GENOME.out.dict,
                file(params.known_snp_vcf, checkIfExists: true),
                file(params.known_snp_vcf_tbi, checkIfExists: true),
                file(params.known_indel_vcf, checkIfExists: true),
                file(params.known_indel_vcf_tbi, checkIfExists: true),
                params.split_reads)
    ch_versions = ch_versions.mix( ALIGN_MARKDUP_BQSR_STATS.out.versions )

    // Somatic variant detection and annotation
    ALIGN_MARKDUP_BQSR_STATS.out.bam.combine(ALIGN_MARKDUP_BQSR_STATS.out.bai, by: 0)
        .branch{ meta, bam, bai ->
            new_meta = meta.clone()
            new_meta.id = meta.pid
            new_meta.pid = ""
            new_meta.tissue = ""
            new_meta.purity = ""
            tumor: bam.name.contains('_T')
                return [new_meta, bam, bai]
           normal: bam.name.contains('_N')
               return [new_meta, bam, bai]
        }
        .set {ch_sample_bams}

    // Estimating Microsatellite instability in tumor samples
    ESTIMATE_MSI ( ch_sample_bams.tumor,
                   Channel.fromPath(params.models, checkIfExists: true).collect() )
    ch_versions = ch_versions.mix(ESTIMATE_MSI.out.versions)
    ESTIMATE_MSI.out.msi
        | MERGE_MSI
        | collectFile (name: 'msi.all.txt', storeDir: "${params.output_dir}/msi")

    // Somatic variant detection and annotation
    ALIGN_MARKDUP_BQSR_STATS.out.cram.combine(ALIGN_MARKDUP_BQSR_STATS.out.crai, by: 0)
        .branch{ meta, cram, crai ->
            new_meta = meta.clone()
            new_meta.id = meta.pid
            new_meta.pid = ""
            new_meta.tissue = ""
            new_meta.purity = ""
            tumor: cram.name.contains('_T')
                return [new_meta, cram, crai]
            normal: cram.name.contains('_N')
                return [new_meta, cram, crai]
        }
        .set {ch_sample_crams}

    ch_sample_crams.tumor.combine(ch_sample_crams.normal, by: 0)
        .map { meta, tumor_cram, tumor_crai, normal_cram, normal_crai ->
            [meta, [tumor_cram, normal_cram], [tumor_crai, normal_crai]] }
        .set { ch_paired_crams }
    bed_files = Channel.fromPath(params.tumor_panel_bed_files, checkIfExists: true)
    ch_paired_crams.combine(bed_files)
        .map { meta, input_crams, input_index_files, intervals ->
            sid = intervals.baseName != "no_intervals" ? new_meta.id + "_" + intervals.baseName : new_meta.id
            new_meta = meta.clone()
            new_meta.sid = intervals.baseName != "no_intervals" ? new_meta.id + "_" + intervals.baseName : new_meta.id
            intervals = intervals.baseName != "no_intervals" ? intervals : []
            [new_meta, input_crams, input_index_files, intervals]
        }
        .set {ch_input_files}
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
    SNV_STRELKA2 (ch_paired_crams,
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