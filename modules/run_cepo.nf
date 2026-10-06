// modules/run_cepo.nf

process RUN_CEPO {
    tag "${h5ad.simpleName}"
    label 'high_mem'
    publishDir "${params.outdir}/metrics/cepo", mode: 'copy'

    container "${projectDir}/environments/py-r-cepo-scdrs.sif"

    input:
    path h5ad

    output:
    path "cepo_results.rds", emit: rds
    path "cepo_norm.csv",    emit: stats
    path "cepo_pvalues.csv", emit: pvalues
    path "*.log",            emit: log

    script:
    """
    Rscript "${projectDir}/bin/run_cepo.R" \\
        --input_h5ad "${h5ad}" \\
        --cell_type_col "${params.cell_type_col}" \\
        --output_rds "cepo_results.rds" \\
        --output_stats "cepo_norm.csv" \\
        --output_pvalues "cepo_pvalues.csv" \\
        --compute_pvalue ${params.cepo_compute_pvalue} \\
        --prefilter_pzero ${params.cepo_prefilter_pzero} \\
        2>&1 | tee run_cepo.log
    """
}
