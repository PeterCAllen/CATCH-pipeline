// modules/mbat_to_seismic_zstat.nf

process MBAT_TO_SEISMIC_ZSTAT {
    tag "${mbat_combined.simpleName}"
    label 'medium_mem'
    publishDir "${params.outdir}/seismic", mode: 'copy'

    container "${projectDir}/environments/py-r-cepo-scdrs.sif"

    input:
    path mbat_combined
    path gene_coords

    output:
    path "${params.gwas_name}.seismic_zstat.tsv", emit: zstat
    path "*.log",                                 emit: log

    script:
    """
    Rscript "${projectDir}/bin/mbat_to_seismic_zstat.R" \\
        "${mbat_combined}" \\
        "${gene_coords}" \\
        "${params.gwas_name}.seismic_zstat.tsv" \\
        2>&1 | tee ${params.gwas_name}.seismic_zstat.log
    """
}
