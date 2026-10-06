// modules/combine_cauchy.nf
// CATCH: Cauchy (ACAT) combination of the two component p-values, per cell type.

process COMBINE_CAUCHY {
    tag "${gwas_name}"
    label 'low_mem'
    publishDir "${params.outdir}/combined", mode: 'copy'

    container "${projectDir}/environments/py-r-cepo-scdrs.sif"

    input:
    tuple path(method_files), val(method_specs)
    val  h5ad_name
    val  gwas_name

    output:
    path "${gwas_name}_catch_combined.tsv", emit: combined_results
    path "${gwas_name}_catch_plot.png",     emit: plot_file
    path "*.log",                           emit: log

    script:
    def spec_args = method_specs.collect { "--method '${it}'" }.join(' \\\n        ')
    """
    Rscript ${projectDir}/bin/cauchy.R \\
        ${spec_args} \\
        --dataset "${h5ad_name}" \\
        --trait "${gwas_name}" \\
        --out "${gwas_name}" \\
        2>&1 | tee ${gwas_name}_cauchy.log
    """
}
