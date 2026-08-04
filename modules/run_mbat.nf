// modules/run_mbat.nf
// Step 1: Run GCTA mBAT per chromosome

process RUN_MBAT {
    tag "chr${chr}"
    label 'medium_mem'
    publishDir "${params.outdir}/mbat", mode: 'copy', pattern: "*_chr*.gene.assoc.mbat"
    
    container "${projectDir}/environments/gcta_v1.94.1.sif"

    input:
    tuple val(chr), path(formatted_gwas), path(mbat_genes), val(plink_prefix), path(plink_files)

    output:
    tuple val(chr), path("*_chr${chr}.gene.assoc.mbat"), emit: mbat_result
    path "*.log", emit: log

    script:
    def gwas_prefix = formatted_gwas.simpleName.replaceAll(~/_formatted$/, '')
    
    """
    echo "========================================" | tee mbat_chr${chr}.log
    echo "Running mBAT for chromosome ${chr}" | tee -a mbat_chr${chr}.log
    echo "========================================" | tee -a mbat_chr${chr}.log
    echo "GWAS: ${formatted_gwas}" | tee -a mbat_chr${chr}.log
    echo "Gene list: ${mbat_genes}" | tee -a mbat_chr${chr}.log
    echo "PLINK prefix: ${plink_prefix}" | tee -a mbat_chr${chr}.log
    echo "Window: ${params.mbat_window_kb} kb" | tee -a mbat_chr${chr}.log
    echo "" | tee -a mbat_chr${chr}.log
    
    gcta64 \\
        --mBAT-combo "${formatted_gwas}" \\
        --bfile "${plink_prefix}.${chr}" \\
        --mBAT-gene-list "${mbat_genes}" \\
        --mBAT-wind ${params.mbat_window_kb} \\
        --chr ${chr} \\
        --out "${gwas_prefix}_chr${chr}" \\
        2>&1 | tee -a mbat_chr${chr}.log
    
    
    N_GENES=\$(tail -n +2 "${gwas_prefix}_chr${chr}.gene.assoc.mbat" | wc -l)
    echo "" | tee -a mbat_chr${chr}.log
    echo "✓ Chr ${chr}: \${N_GENES} genes processed" | tee -a mbat_chr${chr}.log
    """
}
