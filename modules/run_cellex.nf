// modules/run_cellex.nf

process RUN_CELLEX {
    tag "${h5ad.simpleName}"
    label 'high_mem'
    publishDir "${params.outdir}/metrics/cellex", mode: 'copy'

    container "${projectDir}/environments/cellex.sif"

    input:
    path h5ad

    output:
    path "ges.csv",       emit: ges
    path "cellex.*.csv*", emit: all_metrics, optional: true
    path "*.log",         emit: log

    script:
    """
    python3 -u "${projectDir}/bin/run_cellex.py" \\
        --input_h5ad "${h5ad}" \\
        --cell_type_col "${params.cell_type_col}" \\
        --prefix "cellex" \\
        --output_ges "ges.csv" \\
        2>&1 | tee run_cellex.log
    """
}
