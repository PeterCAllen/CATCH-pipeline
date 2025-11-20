// modules/run_scdrs_score.nf
// Step 5: Run scDRS compute-score

process RUN_SCDRS_SCORE {
    tag "${geneset.simpleName}"
    label 'high_mem'
    publishDir "${params.outdir}/scdrs/scores", mode: 'copy'
    
    container "${projectDir}/environments/py-r-cepo-scdrs.sif"

    input:
    tuple path(h5ad), path(geneset), path(cov_file)

    output:
    path "*.full_score.gz", emit: score_file
    path "*.log", emit: log

    script:
    def cov_opt = (cov_file.name != 'NO_FILE') ? "--cov-file ${cov_file}" : ""
    """
    echo "========================================" | tee scdrs_score.log
    echo "Running scDRS compute-score" | tee -a scdrs_score.log
    echo "========================================" | tee -a scdrs_score.log
    echo "H5AD: ${h5ad}" | tee -a scdrs_score.log
    echo "Gene set: ${geneset}" | tee -a scdrs_score.log
    echo "Control gene sets: ${params.scdrs_n_ctrl}" | tee -a scdrs_score.log
    echo "Filter data: ${params.scdrs_filter_data}" | tee -a scdrs_score.log
    echo "Raw count: ${params.scdrs_raw_count}" | tee -a scdrs_score.log
    echo "Covariate file: ${cov_file.name != 'NO_FILE' ? cov_file : 'None'}" | tee -a scdrs_score.log
    echo "" | tee -a scdrs_score.log
    
    export NUMBA_CACHE_DIR=/tmp

    scdrs compute-score \\
        --h5ad-file "${h5ad}" \\
        --h5ad-species human \\
        --gs-file "${geneset}" \\
        --gs-species human \\
        --out-folder . \\
        --n-ctrl ${params.scdrs_n_ctrl} \\
        --flag-filter-data ${params.scdrs_filter_data} \\
        --flag-raw-count ${params.scdrs_raw_count} \\
        ${cov_opt} \\
        2>&1 | tee -a scdrs_score.log
    """
}
