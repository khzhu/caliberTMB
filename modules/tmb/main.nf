//
// TMB_ESTIMATION: Derive TMB estimates from targeted exome capture seqeueuncing data
//

process TMB_CALIBRATION {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(maf_files)
	
    output:
    tuple val(meta), path("*.tsv"), emit: tmb
    tuple val(meta), path("*.maf"), emit: maf
    path "versions.yml"             , emit: versions
	

    //Calculate tumor mutation burden per sample
    script:
    def prefix        = task.ext.prefix ?: "tmb"
    def out_maf       = 'all.maf'
    def input_list = maf_files.collect{ "$it"}.join(' ')

    """
    cat $input_list > $out_maf
    Rscript bin/calculate_tmb.R -m $out_maf -o ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(echo \$(R --version 2>&1)
    END_VERSIONS
    """
}