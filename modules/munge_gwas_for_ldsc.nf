// modules/munge_gwas_for_ldsc.nf

process MUNGE_GWAS_FOR_LDSC {
    tag "${gwas_raw.simpleName}"
    label 'medium_mem'
    publishDir "${params.outdir}/ldsc/munged", mode: 'copy'
    
    container "${projectDir}/environments/ldsc_v1.0.1.sif"

    input:
    path gwas_raw
    val genome_build
    path ref_bim        // BIM file for adding rsIDs

    output:
    path "*.sumstats.gz", emit: gwas_munged
    path "*.log", emit: log

    script:
    def gwas_prefix = gwas_raw.simpleName.replaceAll(~/\.(txt|tsv)$/, '')
    
    """
    echo "========================================" | tee -a munge.log
    echo "STEP 1: Munge GWAS for LDSC" | tee -a munge.log
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
    
    # Step 1b: Munge summary statistics
    echo "" | tee -a munge.log
    echo "Munging summary statistics with LDSC..." | tee -a munge.log
    
    munge_sumstats.py \\
        --sumstats "${gwas_prefix}_with_rsids.txt" \\
        --N-col N \\
        --out "${gwas_prefix}_munged" \\
        2>&1 | tee -a munge.log
    
    echo "" | tee -a munge.log
    echo "✓ Munging complete" | tee -a munge.log
    
    # Verify output
    if [ -f "${gwas_prefix}_munged.sumstats.gz" ]; then
        echo "✓ Output file created: ${gwas_prefix}_munged.sumstats.gz" | tee -a munge.log
        zcat "${gwas_prefix}_munged.sumstats.gz" | head -n 5 | tee -a munge.log
    else
        echo "❌ ERROR: Munged file not created" | tee -a munge.log
        
    fi
    
    mv munge.log ${gwas_prefix}_munge.log
    """
}
