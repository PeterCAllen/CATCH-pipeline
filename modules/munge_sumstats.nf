// modules/munge_sumstats.nf

process MUNGE_SUMSTATS {
    tag "${gwas_formatted.simpleName}"
    label 'medium_mem'
    publishDir "${params.outdir}/conldsc/gwas", mode: 'copy'

    container "${projectDir}/environments/ldsc-py3.sif"

    input:
    path gwas_formatted
    path w_hm3_snplist

    output:
    path "${params.gwas_name}.sumstats.gz", emit: gwas_munged
    path "*.log",                           emit: log

    script:
    """
    munge_sumstats.py \\
        --sumstats "${gwas_formatted}" \\
        --merge-alleles "${w_hm3_snplist}" \\
        --N ${params.gwas_sample_size} \\
        --out "${params.gwas_name}" \\
        2>&1 | tee ${params.gwas_name}.munge_sumstats.stdout.log

    if [ ! -s "${params.gwas_name}.sumstats.gz" ]; then
        echo "ERROR: munge_sumstats.py did not produce ${params.gwas_name}.sumstats.gz"
        exit 1
    fi

    echo "Munged variants: \$(zcat ${params.gwas_name}.sumstats.gz | tail -n +2 | wc -l)"
    """
}
