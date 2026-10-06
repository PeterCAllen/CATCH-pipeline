// modules/cepo_esw_star.nf

process CEPO_ESW_STAR {
    tag "cepo_s"
    label 'medium_mem'
    publishDir "${params.outdir}/metrics/cepo", mode: 'copy'

    container "${projectDir}/environments/py-r-cepo-scdrs.sif"

    input:
    path cepo_stats
    path cepo_pvalues

    output:
    path "cepo_s.csv", emit: stats
    path "*.log",      emit: log

    script:
    """
    Rscript "${projectDir}/bin/esw_star.R" \\
        "${cepo_stats}" \\
        "${cepo_pvalues}" \\
        "cepo_s.csv" \\
        2>&1 | tee cepo_esw_star.log
    """
}
