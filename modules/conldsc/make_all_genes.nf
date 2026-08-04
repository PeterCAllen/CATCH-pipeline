// modules/conldsc/make_all_genes.nf
// CELLECT rule: make_all_genes_background

process MAKE_ALL_GENES_BACKGROUND {
    tag "${specificity_id}"
    label 'low_mem'
    publishDir { "${params.outdir}/conldsc/${specificity_id}/precomputation" }, mode: 'copy'

    container "${projectDir}/environments/cellect-py3.sif"

    input:
    tuple val(specificity_id), path(spec_matrix)

    output:
    tuple val(specificity_id), path("all_genes.${specificity_id}.csv"), emit: all_genes

    script:
    """
    python3 -u "${projectDir}/bin/cellect_make_all_genes.py" \\
        --spec_matrix "${spec_matrix}" \\
        --out "all_genes.${specificity_id}.csv"
    """
}
