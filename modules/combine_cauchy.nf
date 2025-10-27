// modules/combine_cauchy.nf
// Combine p-values from LDSC, MAGMA, and scDRS using Cauchy combination (ACAT)

process COMBINE_CAUCHY {
    tag "${gwas_name}"
    publishDir "${params.outdir}/z_cauchy", mode: 'copy'
    
    container "${projectDir}/environments/py-r-cepo-scdrs.sif"

    input:
    path ldsc_results       // LDSC annotation comparison CSV
    path magma_results      // MAGMA .gsa.out file
    path scdrs_results      // scDRS .scdrs_group file
    val h5ad_name           // Name of single-cell dataset
    val gwas_name           // GWAS trait name/prefix

    output:
    path "${gwas_name}_cauchy_combined.tsv", emit: combined_results
    path "cauchy_combination.log", emit: log
    path "*_cauchy_acato_fdr_plot.png", emit: plot_file

    script:
    """
    echo "========================================" | tee cauchy_combination.log
    echo "Cauchy Combination of P-values" | tee -a cauchy_combination.log
    echo "========================================" | tee -a cauchy_combination.log
    echo "LDSC results  : ${ldsc_results}" | tee -a cauchy_combination.log
    echo "MAGMA results : ${magma_results}" | tee -a cauchy_combination.log
    echo "scDRS results : ${scdrs_results}" | tee -a cauchy_combination.log
    echo "Dataset       : ${h5ad_name}" | tee -a cauchy_combination.log
    echo "GWAS trait    : ${gwas_name}" | tee -a cauchy_combination.log
    echo "" | tee -a cauchy_combination.log

    # Run the Cauchy combination script
    Rscript ${projectDir}/bin/cauchy.R \\
        "${ldsc_results}" \\
        "${magma_results}" \\
        "${scdrs_results}" \\
        "${h5ad_name}" \\
        "${gwas_name}" \\
        2>&1 | tee -a cauchy_combination.log

    echo "" | tee -a cauchy_combination.log
    echo "✓ Cauchy combination complete!" | tee -a cauchy_combination.log
    echo "========================================" | tee -a cauchy_combination.log
    """
}
