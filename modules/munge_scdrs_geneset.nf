// modules/munge_scdrs_geneset.nf
// Step 4: Munge TSV to scDRS gene set format

process MUNGE_SCDRS_GENESET {
    tag "${tsv_file.simpleName}"
    label 'low_mem'
    publishDir "${params.outdir}/scdrs/genesets", mode: 'copy'
    
    container "${projectDir}/environments/py-r-cepo-scdrs.sif"

    input:
    path tsv_file

    output:
    path "*.gs", emit: geneset
    path "*.log", emit: log

    script:
    def gwas_prefix = tsv_file.simpleName
    
    """
    echo "========================================" | tee munge_gs.log
    echo "Creating scDRS gene set file" | tee -a munge_gs.log
    echo "========================================" | tee -a munge_gs.log
    echo "Input: ${tsv_file}" | tee -a munge_gs.log
    echo "" | tee -a munge_gs.log
    
    scdrs munge-gs \\
        --out-file "${gwas_prefix}.gs" \\
        --zscore-file "${tsv_file}" \\
        --weight zscore \\
        --n-max ${params.scdrs_top_genes} \\
        2>&1 | tee -a munge_gs.log
    
    echo "" | tee -a munge_gs.log
    echo "✓ Gene set creation complete" | tee -a munge_gs.log
    echo "  Output: ${gwas_prefix}.gs" | tee -a munge_gs.log
    echo "" | tee -a munge_gs.log
    echo "Gene set preview:" | tee -a munge_gs.log
    head -n 10 "${gwas_prefix}.gs" | tee -a munge_gs.log
    
    mv munge_gs.log ${gwas_prefix}_munge_gs.log
    """
}
