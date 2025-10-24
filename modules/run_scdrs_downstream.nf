// modules/run_scdrs_downstream.nf
// Step 6: Run scDRS downstream analysis (cell type associations)

process RUN_SCDRS_DOWNSTREAM {
    tag "${score_file.simpleName}"
    label 'high_mem'
    publishDir "${params.outdir}/scdrs/results", mode: 'copy'
    
    container "${projectDir}/environments/scDRS_v1.0.4.sif"

    input:
    tuple path(h5ad), path(score_file)

    output:
    path "*group*", emit: group_results
    path "*.log", emit: log

    script:
    """
    echo "========================================" | tee scdrs_downstream.log
    echo "Running scDRS perform-downstream" | tee -a scdrs_downstream.log
    echo "========================================" | tee -a scdrs_downstream.log
    echo "H5AD: ${h5ad}" | tee -a scdrs_downstream.log
    echo "Score file: ${score_file}" | tee -a scdrs_downstream.log
    echo "Cell type column: ${params.cell_type_col}" | tee -a scdrs_downstream.log
    echo "" | tee -a scdrs_downstream.log
    
    scdrs perform-downstream \\
        --h5ad-file "${h5ad}" \\
        --score-file ./"${score_file}" \\
        --out-folder . \\
        --group-analysis "${params.cell_type_col}" \\
        --flag-filter-data ${params.scdrs_filter_data} \\
        --flag-raw-count ${params.scdrs_raw_count} \\
        2>&1 | tee -a scdrs_downstream.log
    """
}
