process GATK4_GETPILEUPSUMMARIES {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(input_bams), path(input_index_files), path(intervals)
    tuple val(meta2), path(fasta)
    tuple val(meta2), path(fai)
    tuple val(meta2), path(dict)
    tuple val(meta3), path(variants), path(variants_tbi)
    val   control_bam

    output:
    tuple val(meta), path('*.pileups.table'), emit: table
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = control_bam ? "${meta.id}.normal" : "${meta.id}.tumor"
    def interval_command = intervals ? "--intervals $intervals" : "--intervals $variants"
    def reference_command = fasta ? "--reference $fasta" : ''
    def bam_file = control_bam? "--input ${input_bams[1]}":"--input ${input_bams[0]}"
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK GetPileupSummaries] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        GetPileupSummaries \\
        $bam_file \\
        --variant $variants \\
        --output ${prefix}.pileups.table \\
        $reference_command \\
        $interval_command \\
        --tmp-dir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pileups.table

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
