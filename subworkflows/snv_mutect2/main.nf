//
// GATK4_MUTECT2: Call variants with Mutect2 (2.1 part of GATK 4)
//

include { GATK4_MUTECT2                                         } from '../../modules/gatk4/mutect2/main'
include { GATK4_GETPILEUPSUMMARIES as GETPILEUPSUMMARIES_NORMAL } from '../../modules/gatk4/getpileupsummaries/main'
include { GATK4_GETPILEUPSUMMARIES as GETPILEUPSUMMARIES_TUMOR  } from '../../modules/gatk4/getpileupsummaries/main'
include { GATK4_CALCULATECONTAMINATION                          } from '../../modules/gatk4/calculatecontamination/main'
include { GATK4_FILTERMUTECTCALLS                               } from '../../modules/gatk4/filtermutectcalls/main'
include { BCFTOOLS_NORM                                         } from '../../modules/bcftools/norm/main'
include { GATK4_MERGEVCFS as MUTECT2_MERGEVCFS                  } from '../../modules/gatk4/mergevcfs/main'
include { VEP as MUTECT2_VEP                                    } from '../../modules/vep/main'
include { VEP as MUTECT2_VEP_JSON                               } from '../../modules/vep/main'
include { VCF2MAF as MUTECT2_VCF2MAF                            } from '../../modules/vcf2maf/main'


workflow SNV_MUTECT2 {

    take:
    ch_input_files                    // channel: [mandatory] [ val(meta), path(ch_input_bams), path(ch_input_bams), path(intervals) ]
    ch_fasta                          // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_fai                            // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_dict                           // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_germline_resource
    ch_pileup_variants
    ch_cosmic_vcf
    vep_cache

    main:
    ch_versions         = Channel.empty()

    GATK4_MUTECT2 ( ch_input_files,
                    ch_fasta, ch_fai, ch_dict,
                    ch_germline_resource)
    ch_versions = ch_versions.mix(GATK4_MUTECT2.out.versions)

    GETPILEUPSUMMARIES_TUMOR  ( ch_input_files,
                                ch_fasta, ch_fai, ch_dict,
                                ch_pileup_variants,
                                false)

    GETPILEUPSUMMARIES_NORMAL ( ch_input_files,
                                ch_fasta, ch_fai, ch_dict,
                                ch_pileup_variants,
                                true)

    GATK4_CALCULATECONTAMINATION ( GETPILEUPSUMMARIES_TUMOR.out.table,
                                    GETPILEUPSUMMARIES_NORMAL.out.table)
    GATK4_FILTERMUTECTCALLS ( GATK4_MUTECT2.out.vcf.combine(GATK4_MUTECT2.out.tbi,
                                by:0).combine(GATK4_MUTECT2.out.stats, by:0),
                                GATK4_CALCULATECONTAMINATION.out.contamination,
                                ch_fasta, ch_fai, ch_dict)
    
    BCFTOOLS_NORM ( GATK4_FILTERMUTECTCALLS.out.vcf.combine(GATK4_FILTERMUTECTCALLS.out.tbi, by:0), ch_fasta )
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions)

    ch_vcf_files = BCFTOOLS_NORM
                        .out.vcf
                        .map { meta, norm_vcfs ->
                            def new_meta = [id: meta.pid]
                            [new_meta, norm_vcfs] }
                        .groupTuple()
    MUTECT2_MERGEVCFS ( ch_vcf_files, ch_dict )

    MUTECT2_VEP ( MUTECT2_MERGEVCFS.out.vcf.combine(MUTECT2_MERGEVCFS.out.tbi, by:0),
                    ch_fasta,
                    ch_germline_resource,
                    ch_cosmic_vcf,
                    vep_cache)
    ch_versions = ch_versions.mix(MUTECT2_VEP.out.versions)

    MUTECT2_VEP_JSON ( MUTECT2_MERGEVCFS.out.vcf.combine(MUTECT2_MERGEVCFS.out.tbi, by:0),
                    ch_fasta,
                    ch_germline_resource,
                    ch_cosmic_vcf,
                    vep_cache)

    MUTECT2_VCF2MAF ( MUTECT2_VEP.out.vcf, ch_fasta, vep_cache)
    ch_versions = ch_versions.mix(MUTECT2_VCF2MAF.out.versions)

    emit:
    vcf      = MUTECT2_MERGEVCFS.out.vcf    // channel: [ val(meta), path("*.vcf.gz") ]
    tbi      = MUTECT2_MERGEVCFS.out.tbi    // channel: [ val(meta), path("*.tbi") ]
    vep      = MUTECT2_VEP.out.vcf          // channel: [ val(meta), path("*.vep.vcf") ]
    html     = MUTECT2_VEP.out.html         // channel: [ val(meta), path("*.html")]
    json     = MUTECT2_VEP_JSON.out.json    // channel: [ val(meta), path("*.vep.json") ]
    maf      = MUTECT2_VCF2MAF.out.maf      // channel: [ val(meta), path("*.maf") ]
    stats    = GATK4_MUTECT2.out.stats      // channel: [ val(meta), path("*.stats") ]
    f1r2     = GATK4_MUTECT2.out.f1r2
    versions = ch_versions                  // channel: [ path(versions.yml) ]
}