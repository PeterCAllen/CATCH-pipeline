// modules/run_seismic.nf

process RUN_SEISMIC {
    tag "${params.gwas_name}"
    label 'high_mem'
    publishDir "${params.outdir}/seismic", mode: 'copy'

    container "${projectDir}/environments/seismic.sif"

    input:
    tuple path(h5ad), path(zstat)

    output:
    path "${params.gwas_name}_seismic.tsv",            emit: results
    path "${params.gwas_name}_seismic_specificity.csv", emit: specificity
    path "*.log",                                       emit: log

    script:
    """
    Rscript "${projectDir}/bin/run_seismic.R" \\
        --input_h5ad "${h5ad}" \\
        --cell_type_col "${params.cell_type_col}" \\
        --zstat_file "${zstat}" \\
        --output_tsv "${params.gwas_name}_seismic.tsv" \\
        --output_sscore "${params.gwas_name}_seismic_specificity.csv" \\
        --assay_name "${params.seismic_assay}" \\
        --lognorm ${params.seismic_lognorm} \\
        2>&1 | tee ${params.gwas_name}_seismic.log
    """
}
