process GATK4_MUTECT2 {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(input_bams), path(input_index_files), path(intervals)
    tuple val(meta2), path(fasta)
    tuple val(meta2), path(fai)
    tuple val(meta2), path(dict)
    tuple val(meta3), path(germline_resource), path(germline_resource_tbi)

    output:
    tuple val(meta), path("*.vcf.gz")     , emit: vcf
    tuple val(meta), path("*.tbi")        , emit: tbi
    tuple val(meta), path("*.stats")      , emit: stats
    tuple val(meta), path("*.f1r2.tar.gz"), optional:true, emit: f1r2
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.mutect2"
    def inputs = input_bams.collect{ "--input $it"}.join(" ")
    def interval_command = intervals ? "--intervals ${intervals}" : ""
    def tumor_sample = "--tumor-sample ${meta.pid}_T"
    def normal_sample = "--normal-sample ${meta.pid}_N"
    def gr_command = germline_resource ? "--germline-resource $germline_resource" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK Mutect2] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        Mutect2 \\
        $args \\
        $inputs \\
        --native-pair-hmm-threads ${task.cpus} \\
        --output ${prefix}.vcf.gz \\
        $tumor_sample \\
        $normal_sample \\
        --reference $fasta \\
        $gr_command \\
        --f1r2-tar-gz ${prefix}.f1r2.tar.gz \\
        $interval_command \\
        --tmp-dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    touch ${prefix}.vcf.gz.stats
    touch ${prefix}.f1r2.tar.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
