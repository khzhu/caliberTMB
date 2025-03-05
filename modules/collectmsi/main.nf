//
// MERGE_MSI: Gather MSI score for each sample
//

process MERGE_MSI {
    tag "${params.batch_id}"
    label 'process_small'

    input:
    path(msi_text)
    output:
    path('msi_all.txt')

    script:
    """
    paste -d'\t' <(basename ${msi_text} |cut -d"." -f1 |cut -d"_" -f1) <(tail -1 ${msi_text}) > msi_all.txt

    """
}