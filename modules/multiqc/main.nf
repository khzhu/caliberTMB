process MULTIQC {
    tag "${params.batch_id}"
    label 'process_single'

    container = "${params.container_dir}/multiqc_v1.28.simg"
    cpus = 16
    memory = 64.GB
    ext.args = '-m mosdepth -m picard -m samtools -m fastqc'
    publishDir (
            path: { "${params.output_dir}/cohort/" },
            mode: params.publish_dir_mode,
            pattern: "*{html,yml,_data}",
            saveAs: { "multiqc/${it}" }
    )

    input:
    path bam_files

    output:
    path "*multiqc_report.html", emit: report
    path "*multiqc_data"       , emit: data
    path "*_plots"             , optional:true, emit: plots
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ? "--filename ${task.ext.prefix}.html" : "--filename ${params.batch_id}.html"
    def subdirs = new File("${params.output_dir}").listFiles().findAll {
                    it.isDirectory() && (it.name.contains('_T') || it.name.contains('_N'))}
                    .collect { it.absolutePath + "/qc" }.join(" ")

    """
    multiqc \\
        --force \\
        $args \\
        $subdirs \\
        -o .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
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