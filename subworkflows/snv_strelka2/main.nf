//
// GATK4_MUTECT2: Call variants with Mutect2 (2.1 part of GATK 4)
//

include { MANTA_SOMATIC                        } from '../../modules/manta/main'
include { STRELKA_SOMATIC                      } from '../../modules/strelka/main'
include { BCFTOOLS_NORM as BCFTOOLS_NORM_SNV   } from '../../modules/bcftools/norm/main'
include { BCFTOOLS_NORM as BCFTOOLS_NORM_INDEL } from '../../modules/bcftools/norm/main'
include { GATK4_MERGEVCFS                      } from '../../modules/gatk4/mergevcfs/main'

workflow SNV_STRELKA2 {

    take:
    ch_input_bams   // channel: [mandatory] [ val(meta), path(bam_tumor), path(bam_normal) ]
    ch_fasta        // channel: [mandatory] [ val(meta2), path(fasta) ]
    ch_fai          // channel: [mandatory] [ val(meta2), path(fai) ]
    ch_dict         // channel: [mandatory] [ val(meta2), path(dict) ]
    ch_target_bed   // channel: [mandatory] [ val(meta3), path(bed), path(bed_tbi) ]

    main:
    ch_versions         = Channel.empty()

    MANTA_SOMATIC ( ch_input_bams, 
                    ch_fasta, ch_fai,
                    ch_target_bed)
    ch_versions = ch_versions.mix(MANTA_SOMATIC.out.versions)

    STRELKA_SOMATIC  ( ch_input_bams,
                        MANTA_SOMATIC.out.candidate_small_indels_vcf.combine(
                        MANTA_SOMATIC.out.candidate_small_indels_vcf_tbi, by:0),
                        ch_fasta, ch_fai, ch_target_bed)
    ch_versions = ch_versions.mix(STRELKA_SOMATIC.out.versions)

    BCFTOOLS_NORM_SNV  ( STRELKA_SOMATIC.out.vcf_snvs.combine(
                        STRELKA_SOMATIC.out.vcf_snvs_tbi, by:0),
                        ch_fasta)
    ch_versions = ch_versions.mix(BCFTOOLS_NORM_SNV.out.versions)

    BCFTOOLS_NORM_INDEL ( STRELKA_SOMATIC.out.vcf_indels.combine(
                        STRELKA_SOMATIC.out.vcf_indels_tbi, by:0),
                        ch_fasta)

    GATK4_MERGEVCFS ( BCFTOOLS_NORM_SNV.out.vcf.join(BCFTOOLS_NORM_INDEL.out.vcf), ch_dict )
    ch_versions = ch_versions.mix(GATK4_MERGEVCFS.out.versions)

    emit:
    vcf      = GATK4_MERGEVCFS.out.vcf       // channel: [ val(meta), path("*.vcf.gz") ]
    tbi      = GATK4_MERGEVCFS.out.tbi       // channel: [ val(meta), path("*.tbi") ]
    versions = ch_versions                   // channel: [ path(versions.yml) ]
}
