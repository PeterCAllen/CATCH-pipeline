// modules/format_gwas_for_scdrs.nf
// Step 0: Format GWAS summary statistics for mBAT

process FORMAT_GWAS_FOR_SCDRS {
    tag "${gwas_raw.simpleName}"
    label 'low_mem'
    publishDir "${params.outdir}/scdrs/gwas", mode: 'copy', pattern: "*_formatted.txt"
    
    container "${projectDir}/environments/py-r-cepo-scdrs.sif"

    input:
    path gwas_raw
    val genome_build
    path ref_bim

    output:
    path "*_formatted.txt", emit: formatted_gwas
    path "*_with_rsids.txt", emit: gwas_with_rsids
    path "*.log", emit: log

    script:
    def gwas_prefix = gwas_raw.simpleName.replaceAll(~/\\.txt$/, '').replaceAll(~/\\.gz$/, '')
    
    """
    echo "========================================" | tee format_gwas.log
    echo "STEP 0: Format GWAS for mBAT" | tee -a format_gwas.log
    echo "========================================" | tee -a format_gwas.log
    echo "Input: ${gwas_raw}" | tee -a format_gwas.log
    echo "Genome build: ${genome_build}" | tee -a format_gwas.log
    echo "Reference BIM: ${ref_bim}" | tee -a format_gwas.log
    echo "Prefix: ${gwas_prefix}" | tee -a format_gwas.log
    echo "" | tee -a format_gwas.log
    
    # Step 0a: Add rsIDs to GWAS file using reference panel
    echo "Adding rsIDs from reference panel..." | tee -a format_gwas.log
    bash ${projectDir}/bin/format_gwas_for_ldsc.sh \\
        "${gwas_raw}" \\
        "${ref_bim}" \\
        "${gwas_prefix}" \\
        2>&1 | tee -a format_gwas.log
    
    # Step 0b: Format for mBAT (SNP, A1, A2, freq, b, se, p, N)
    echo "" | tee -a format_gwas.log
    echo "Formatting for mBAT..." | tee -a format_gwas.log
    
    Rscript - "${gwas_prefix}_with_rsids.txt" "${gwas_prefix}_formatted.txt" <<'RSCRIPT'
    library(tidyverse)
    
    args <- commandArgs(trailingOnly = TRUE)
    gwas_file <- args[1]
    out_file <- args[2]
    
    cat(sprintf("Reading GWAS file with rsIDs: %s\\n", gwas_file))
    
    # Read the file with rsIDs (space-delimited from format_gwas_for_ldsc.sh)
    gwas <- read_delim(gwas_file, delim = " ", show_col_types = FALSE)
    
    cat(sprintf("  Columns found: %s\\n", paste(colnames(gwas), collapse=", ")))
    
    # Check if we have required columns
    # Expected columns: SNP, BP, SNP_POS, A1, A2, b, se, P, freq, N
    if (!all(c("SNP", "A1", "A2", "b", "se", "P", "N") %in% colnames(gwas))) {
        # Try alternative column names
        if ("p" %in% colnames(gwas)) gwas <- gwas %>% rename(P = p)
        if ("beta" %in% colnames(gwas)) gwas <- gwas %>% rename(b = beta)
    }
    
    # Add freq if missing (use 0.5 as placeholder if not available)
    if (!"freq" %in% colnames(gwas)) {
        gwas <- gwas %>% mutate(freq = 0.5)
    }
    
    # Format for mBAT: SNP, A1, A2, freq, b, se, p, N
    gwas_formatted <- gwas %>%
      select(SNP, A1, A2, freq, b, se, p = P, N)
    
    write_delim(gwas_formatted, out_file, delim = "\\t")
    
    cat(sprintf("✓ Formatted %d SNPs\\n", nrow(gwas_formatted)))
    cat(sprintf("  Output: %s\\n", out_file))
    cat(sprintf("\\nFirst 5 SNPs:\\n"))
    print(head(gwas_formatted, 5))
    RSCRIPT
    
    if [ ! -f "${gwas_prefix}_formatted.txt" ]; then
        echo "❌ ERROR: Failed to format GWAS file" | tee -a format_gwas.log
        
    fi
    
    echo "" | tee -a format_gwas.log
    echo "✓ Formatting complete" | tee -a format_gwas.log
    
    mv format_gwas.log ${gwas_prefix}_format.log
    """
}
