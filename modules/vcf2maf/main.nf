process VCF2MAF {

    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(input_vep_vcf) // Use an uncompressed VCF file!
    tuple val(meta2), path(fasta)    // Required
    path vep_cache                   // Required for VEP running. A default of /.vep is supplied.

    output:
    tuple val(meta), path("*.maf"), emit: maf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args   ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def tumor_id = task.tumorID         ?: "${meta.id}_T"
    def normal_id = task.normalID       ?: "${meta.id}_N"
    // If VEP is present, it will find it and add it to commands.
    // If VEP is not present they will be blank
    def VERSION = '1.6.22' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    vcf2maf.pl \\
        $args \\
        --ref-fasta $fasta \\
        --input-vcf $input_vep_vcf \\
        --tumor-id ${meta.id}_T  \\
        --normal-id ${meta.id}_N \\
        --vcf-tumor-id $tumor_id \\
        --vcf-normal-id $normal_id \\
        --output-maf ${prefix}.maf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcf2maf: $VERSION
    END_VERSIONS
    """
}
