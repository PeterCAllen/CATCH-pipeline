// modules/sanitize_specificity.nf

process SANITIZE_SPECIFICITY {
    tag "${specificity_id}"
    label 'low_mem'
    publishDir "${params.outdir}/metrics/sanitized", mode: 'copy'

    container "${projectDir}/environments/py-r-cepo-scdrs.sif"

    input:
    tuple val(specificity_id), path(matrix), val(input_sep)

    output:
    tuple val(specificity_id), path("${specificity_id}.specificity.csv"), emit: matrix
    tuple val(specificity_id), path("${specificity_id}.annotations.txt"), emit: annotations
    path "*.log",                                                          emit: log

    script:
    def top_pct = params.conldsc_top_pct ? "--top_pct ${params.conldsc_top_pct}" : ''
    """
    python3 -u "${projectDir}/bin/sanitize_specificity_matrix.py" \\
        --input "${matrix}" \\
        --output "${specificity_id}.specificity.csv" \\
        --output_annotations "${specificity_id}.annotations.txt" \\
        --input_sep ${input_sep} \\
        ${top_pct} \\
        2>&1 | tee ${specificity_id}.sanitize.log
    """
}
