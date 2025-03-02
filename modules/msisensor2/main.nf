process MSISENSOR2_MSI {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bam_index)
    path(models, stageAs: "models/*")

    output:
    tuple val(meta), path("*msi.txt")      , emit: msi
    tuple val(meta), path("*dis.txt")      , emit: distribution
    tuple val(meta), path("*somatic.txt")  , emit: somatic
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def model_cmd     = models     ? "-M models/*"     : ""
    def tumor_bam_cmd = tumor_bam  ? "-t $tumor_bam"  : ""
    """
    /usr/local/bin/msisensor2 msi ${task.ext.args} \\
        -b ${task.cpus} \\
        $model_cmd \\
        $tumor_bam_cmd \\
        -o $prefix 2>&1

    mv ${prefix} ${prefix}_msi.txt
    mv ${prefix}_dis ${prefix}_dis.txt
    mv ${prefix}_somatic ${prefix}_somatic.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        msisensor2: \$(echo \$(/usr/local/bin/msisensor2 2> >(grep Version) | sed 's/Version: v//g'))
    END_VERSIONS
    """
}