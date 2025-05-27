process MULTIQC {
    tag "${params.batch_id}"
    label 'process_single'

    input:
    path  multiqc_path

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ? "--filename ${task.ext.prefix}.html" : "${params.batch_id}.html"
    def subdirs = new File(multiqc_path).listFiles().findAll {
                    it.isDirectory() && (it.name.contains('_T') || it.name.contains('_N'))}
                    .collect { it.absolutePath + "/qc" }

    """
    /usr/local/bin/multiqc \\
        --force \\
        $args \\
        $prefix \\
        $subdirs.join(' ') \\
        -o .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( /usr/local/bin/multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """

    stub:
    """
    mkdir multiqc_data
    mkdir multiqc_plots
    touch multiqc_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}