// modules/munge_gwas_for_ldsc.nf

process CONVERT_GWAS_FOR_LDSC {
    tag "${gwas_raw.simpleName}"
    label 'medium_mem'
    publishDir "${params.outdir}/conldsc/gwas", mode: 'copy'

    container "${projectDir}/environments/py-r-cepo-scdrs.sif"

    input:
    path gwas_raw
    val genome_build
    path ref_bim

    output:
    path "${gwas_raw.simpleName}_for_ldsc.txt", emit: formatted_gwas
    path "*.log",                               emit: log

    script:
    def gwas_prefix = gwas_raw.simpleName
    """
    bash ${projectDir}/bin/format_gwas_for_ldsc.sh \\
        "${gwas_raw}" \\
        "${ref_bim}" \\
        "${gwas_prefix}" \\
        2>&1 | tee ${gwas_prefix}.rsid.log

    Rscript ${projectDir}/bin/munge_gwas_for_ldsc.R \\
        "${gwas_prefix}_with_rsids.txt" \\
        "${gwas_prefix}_for_ldsc.txt" \\
        "${params.gwas_sample_size}" \\
        2>&1 | tee ${gwas_prefix}.harmonize.log
    """
}
