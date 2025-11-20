// modules/munge_gwas_for_ldsc.nf
// Step 1: Munge GWAS summary statistics for LDSC

process CONVERT_GWAS_FOR_LDSC {
    tag "${gwas_raw.simpleName}"
    label 'medium_mem'
    
    container "${projectDir}/environments/py-r-cepo-scdrs.sif"

    input:
    path gwas_raw
    val genome_build
    path ref_bim

    output:
    path "*_for_ldsc.txt"

    script:
    def gwas_prefix = gwas_raw.simpleName.replaceAll(~/\\.txt$/, '').replaceAll(~/\\.gz$/, '')
    
    """
    echo "========================================" | tee -a munge.log
    echo "STEP 1a: Format GWAS for LDSC" | tee -a munge.log
    echo "========================================" | tee -a munge.log
    echo "Genome build: ${genome_build}" | tee -a munge.log
    echo "Reference BIM: ${ref_bim}" | tee -a munge.log
    echo "" | tee -a munge.log
    
    # Step 1a: Add rsIDs to GWAS file
    echo "Adding rsIDs from reference panel..." | tee -a munge.log
    bash ${projectDir}/bin/format_gwas_for_ldsc.sh \\
        "${gwas_raw}" \\
        "${ref_bim}" \\
        "${gwas_prefix}"
    
    # Step 1b: Convert Z-scores to P-values if needed
    echo "" | tee -a munge.log
    echo "Checking for Z-scores and converting to P-values if needed..." | tee -a munge.log
    
    Rscript ${projectDir}/bin/munge_gwas_for_ldsc.R "${gwas_prefix}_with_rsids.txt" "${gwas_prefix}_for_ldsc.txt"
    """
}
