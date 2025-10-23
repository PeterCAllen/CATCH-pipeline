// modules/local/run_cepo.nf

process RUN_CEPO {
    tag "$h5ad_processing"    
    container "${projectDir}/environments/py-r-cepo-scdrs.sif"
    publishDir "${params.outdir}/cepo", mode: 'copy'

    input:
    path h5ad_processed
    path cepo_coords

    output:
    path "cepo_results.rds", emit: rds
    path "cepo_stats.tsv", emit: stats_tsv

    script:
    """
    Rscript "${projectDir}/bin/run_cepo.R" \\
        --input_h5ad "${h5ad_processed}" \\
        --cell_type_col "${params.cell_type_col}" \\
        --output_rds cepo_results.rds \\
        --output_tsv cepo_stats.tsv \\
        --gene_coords "${cepo_coords}"
    """
}
