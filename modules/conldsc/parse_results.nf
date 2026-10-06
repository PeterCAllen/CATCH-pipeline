// modules/conldsc/parse_results.nf
// CELLECT rule: parse_results (prioritization block)
//
// The work-dir filename is prefixed with the specificity id so that several
// specificity matrices can be staged side by side; the published copy keeps
// CELLECT's plain prioritization.csv name for direct diffing.

process PARSE_LDSC_RESULTS {
    tag "${specificity_id}"
    label 'low_mem'
    publishDir { "${params.outdir}/conldsc/${specificity_id}/results" },
               mode: 'copy',
               saveAs: { 'prioritization.csv' }

    container "${projectDir}/environments/cellect-py3.sif"

    input:
    tuple val(specificity_id), path(cell_type_results)

    output:
    tuple val(specificity_id), path("${specificity_id}_prioritization.csv"), emit: prioritization

    script:
    """
    python3 -u "${projectDir}/bin/cellect_parse_results.py" \\
        --cell_type_results "${cell_type_results}" \\
        --gwas "${params.gwas_name}" \\
        --out "${specificity_id}_prioritization.csv"
    """
}
