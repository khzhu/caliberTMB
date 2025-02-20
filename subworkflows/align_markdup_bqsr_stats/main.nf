//
// Alignment to genome with BWA and sort reads by SAMBAMBA
//

include { SEQKIT_SPLIT2              } from '../../modules/seqkit/main'
include { BWA_MEM as BWA_MEM1        } from '../../modules/bwa/mem/main'
include { BWA_MEM as BWA_MEM2        } from '../../modules/bwa/mem/main'
include { SAMBAMBA_MERGE             } from '../../modules/sambamba/merge/main'
include { SAMTOOLS_SORT as SAM_SORT1 } from '../../modules/samtools/sort/main'
include { SAMTOOLS_SORT as SAM_SORT2 } from '../../modules/samtools/sort/main'
include { SAMTOOLS_SORT as SAM_SORT  } from '../../modules/samtools/sort/main'
include { SAMBAMBA_MARKDUP           } from '../../modules/sambamba/markdup/main'
include { SAMBAMBA_FLAGSTAT          } from '../../modules/sambamba/flagstat/main'
include { GATK4_BASERECALIBRATOR     } from '../../modules/gatk4/baserecalibrator/main'
include { GATK4_APPLYBQSR            } from '../../modules/gatk4/applybqsr/main'
include { SAMTOOLS_STATS             } from '../../modules/samtools/stats/main'
include { SAMTOOLS_CONVERT           } from '../../modules/samtools/convert/main'


workflow ALIGN_MARKDUP_BQSR_STATS {
    take:
    reads                  // channel (mandatory): [ val(meta), [ path(reads) ] ]
    bwa_index              // channel (mandatory): [ val(meta2), [ path(index) ] ]
    fasta                  // channel (optional) : [ val(meta2), [ path(fasta) ] ]
    val_sort_bam           // boolean (mandatory): true or false
    intervals              // channel: [ path(intervals) ]
    fai                    // channel: [ val(meta2), path(fai) ]
    dict                   // channel: [ val(meta2), path(dict) ]
    snp_known_sites        // channel: [ path(known_sites) ]
    snp_known_sites_tbi    // channel: [ path(known_sites_tbi) ]
    indel_known_sites      // channel: [ path(known_sites) ]
    indel_known_sites_tbi  // channel: [ path(known_sites_tbi)
    split_reads            // boolean (mandatory): true or false

    main:
    ch_versions = Channel.empty()

    //
    // Map reads with BWA
    //
    if ( split_reads ) {
        //Split trimmed reads into 2 parts
        SEQKIT_SPLIT2 ( reads )
        ch_versions = ch_versions.mix(SEQKIT_SPLIT2.out.versions)
        BWA_MEM1 ( SEQKIT_SPLIT2.out.reads_part1, bwa_index, [[id:'genome'],fasta], val_sort_bam )
        BWA_MEM2 ( SEQKIT_SPLIT2.out.reads_part2, bwa_index, [[id:'genome'],fasta], val_sort_bam )
        ch_versions = ch_versions.mix(BWA_MEM1.out.versions)
        SAM_SORT1 (BWA_MEM1.out.bam, [[id:'genome'],fasta])
        SAM_SORT2 (BWA_MEM2.out.bam, [[id:'genome'],fasta])

        SAMBAMBA_MERGE ( SAM_SORT1.out.bam.combine(SAM_SORT2.out.bam, by:0)
                            .map { it -> tuple(it[0],[it[1],it[2]])} )
        // Sort bam with samtools
        SAM_SORT ( SAMBAMBA_MERGE.out.bam , [[id:'genome'],fasta] )
        ch_versions = ch_versions.mix(SAM_SORT.out.versions)
    } else {
        BWA_MEM1 ( reads, bwa_index, [[id:'genome'],fasta], val_sort_bam )
        ch_versions = ch_versions.mix(BWA_MEM1.out.versions)
        // Sort bam with samtools
        SAM_SORT ( BWA_MEM1.out.bam , [[id:'genome'],fasta] )
        ch_versions = ch_versions.mix(SAM_SORT.out.versions)
    }

    //
    // Run sambamba deduplicate and flagstat
    //
    SAMBAMBA_MARKDUP ( SAM_SORT.out.bam )
    ch_versions = ch_versions.mix(SAMBAMBA_MARKDUP.out.versions)

    //
    // Generate sambamba flagstat
    //
    SAMBAMBA_FLAGSTAT ( SAMBAMBA_MARKDUP.out.bam )

    //
    // Apply for gatk4 BQSR
    //
    GATK4_BASERECALIBRATOR ( SAMBAMBA_MARKDUP.out.bam,
                            intervals, [[id:'genome'],fasta], fai, dict,
                            snp_known_sites, snp_known_sites_tbi,
                            indel_known_sites, indel_known_sites_tbi)
    GATK4_APPLYBQSR ( SAMBAMBA_MARKDUP.out.bam.combine(
                        GATK4_BASERECALIBRATOR.out.table, by: 0),
                        intervals, [[id:'genome'],fasta], fai, dict)
    ch_versions = ch_versions.mix(GATK4_APPLYBQSR.out.versions)

    //
    // Generate samtool stats
    //
    SAMTOOLS_STATS ( GATK4_APPLYBQSR.out.bam, [[id:'genome'],fasta])

    //
    // Convert to CRAM
    //
    SAMTOOLS_CONVERT ( GATK4_APPLYBQSR.out.bam, [[id:'genome'],fasta])

    emit:
    cram     = SAMTOOLS_CONVERT.out.bam       // channel: [ val(meta), path(cram) ]
    crai     = SAMTOOLS_CONVERT.out.crai      // channel: [ val(meta), path(crai) ]
    flagstat = SAMBAMBA_FLAGSTAT.out.stats    // channel: [ val(meta), path(flagstat) ]
    stats    = SAMTOOLS_STATS.out.stats       // channel: [ val(meta), path(stats) ]
    versions = ch_versions                    // channel: [ path(versions.yml) ]
}