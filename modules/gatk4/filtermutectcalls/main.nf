process GATK4_FILTERMUTECTCALLS {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta),  path(input_vcf), path(input_tbi), path(input_stats)
    tuple val(meta),  path(table)
    tuple val(meta2), path(fasta)
    tuple val(meta2), path(fai)
    tuple val(meta2), path(dict)

    output:
    tuple val(meta), path("*.vcf.gz")            , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi")        , emit: tbi
    tuple val(meta), path("*.filteringStats.tsv"), emit: stats
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.filtered"
    def table_command = table ? table.collect{"--contamination-table $it"}.join(' ') : ''

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK FilterMutectCalls] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        FilterMutectCalls \\
        --variant $input_vcf \\
        --output ${prefix}.vcf.gz \\
        --reference $fasta \\
        $table_command \\
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
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    touch ${prefix}.vcf.gz.filteringStats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
