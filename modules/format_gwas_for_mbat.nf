// modules/format_gwas_for_mbat.nf
// Step 0: Format GWAS summary statistics for mBAT

process FORMAT_GWAS_FOR_MBAT {
    tag "${gwas_raw.simpleName}"
    label 'low_mem'
    publishDir "${params.outdir}/mbat/gwas", mode: 'copy', pattern: "*_formatted.txt"
    
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
    
    Rscript - "${gwas_prefix}_with_rsids.txt" "${gwas_prefix}_formatted.txt" <<'EOF'
    library(tidyverse)
    
    args <- commandArgs(trailingOnly = TRUE)
    gwas_file <- args[1]
    out_file <- args[2]
    
    cat(sprintf("Reading GWAS file with rsIDs: %s\\n", gwas_file))
    
    # Simple direct approach - try each delimiter and pick best one
    gwas <- NULL
    delim_name <- ""
    max_ncol <- 0
    
    # Try whitespace
    tryCatch({
        test <- read_table(gwas_file, n_max = 2, show_col_types = FALSE)
        if (ncol(test) > max_ncol) {
            max_ncol <- ncol(test)
            delim_name <- "WHITESPACE"
            gwas <- read_table(gwas_file, show_col_types = FALSE)
        }
    }, error = function(e) {})
    
    # Try tab
    tryCatch({
        test <- read_delim(gwas_file, delim = "\\t", n_max = 2, show_col_types = FALSE)
        if (ncol(test) > max_ncol) {
            max_ncol <- ncol(test)
            delim_name <- "TAB"
            gwas <- read_delim(gwas_file, delim = "\\t", show_col_types = FALSE)
        }
    }, error = function(e) {})
    
    # Try comma
    tryCatch({
        test <- read_delim(gwas_file, delim = ",", n_max = 2, show_col_types = FALSE)
        if (ncol(test) > max_ncol) {
            max_ncol <- ncol(test)
            delim_name <- "COMMA"
            gwas <- read_delim(gwas_file, delim = ",", show_col_types = FALSE)
        }
    }, error = function(e) {})
    
    if (is.null(gwas)) {
        stop("Could not parse file with any delimiter")
    }
    
    cat(sprintf("  Detected delimiter: %s (%d columns)\\n", delim_name, ncol(gwas)))
    cat(sprintf("  Columns found: %s\\n", paste(colnames(gwas), collapse=", ")))
    cat(sprintf("  Total rows: %d\\n", nrow(gwas)))
    
    # Check if we have required columns
    required_cols <- c("SNP", "A1", "A2", "N")
    missing_required <- setdiff(required_cols, colnames(gwas))
    
    if (length(missing_required) > 0) {
        stop(sprintf("ERROR: Missing required columns: %s", paste(missing_required, collapse=", ")))
    }
    
    # Handle different column name variations
    if ("p" %in% colnames(gwas) && !"P" %in% colnames(gwas)) {
        gwas <- gwas %>% rename(P = p)
    }
    if ("beta" %in% colnames(gwas) && !"b" %in% colnames(gwas)) {
        gwas <- gwas %>% rename(b = beta)
    }
    if ("BETA" %in% colnames(gwas) && !"b" %in% colnames(gwas)) {
        gwas <- gwas %>% rename(b = BETA)
    }
    
    # Check if we have effect size columns or Z-scores
    has_beta <- "b" %in% colnames(gwas)
    has_se <- "se" %in% colnames(gwas) || "SE" %in% colnames(gwas)
    has_pval <- "P" %in% colnames(gwas)
    has_z <- "Z" %in% colnames(gwas)
    
    if ("SE" %in% colnames(gwas)) gwas <- gwas %>% rename(se = SE)
    
    # Convert Z-scores if that's what we have
    if (has_z && !has_pval) {
        cat("  Converting Z-scores to p-values\\n")
        gwas <- gwas %>% mutate(P = 2 * pnorm(-abs(Z)))
    }
    
    if (has_z && !has_beta) {
        cat("  WARNING: No effect size (beta) column found. Using Z-scores as proxy.\\n")
        gwas <- gwas %>% mutate(b = Z)
    }
    
    if (has_z && !has_se) {
        cat("  WARNING: No standard error column found. Estimating from sample size.\\n")
        gwas <- gwas %>% mutate(se = abs(Z / qnorm(P/2)))
    }
    
    # Verify we now have all needed columns
    if (!all(c("b", "se", "P") %in% colnames(gwas))) {
        missing <- setdiff(c("b", "se", "P"), colnames(gwas))
        stop(sprintf("ERROR: Missing columns needed for mBAT: %s. Cannot proceed without effect size, standard error, and p-value.", paste(missing, collapse=", ")))
    }
    
    # Add freq if missing (use 0.5 as placeholder if not available)
    if (!"freq" %in% colnames(gwas)) {
        cat("  Adding freq column (using 0.5 as placeholder)\\n")
        gwas <- gwas %>% mutate(freq = 0.5)
    }
    
    # Format for mBAT: SNP, A1, A2, freq, b, se, p, N
    gwas_formatted <- gwas %>%
      select(SNP, A1, A2, freq, b, se, p = P, N) %>%
      filter(!is.na(SNP), !is.na(N), N > 0) %>%
      filter(grepl("^rs", SNP))
    
    # Check if we have any data left
    if (nrow(gwas_formatted) == 0) {
        stop("ERROR: No SNPs with valid rsIDs found after filtering. Check that SNP IDs were successfully converted.")
    }
    
    write_delim(gwas_formatted, out_file, delim = "\\t")
    
    cat(sprintf("✓ Formatted %d SNPs (%.1f%% of input retained)\\n", 
                nrow(gwas_formatted), 
                100 * nrow(gwas_formatted) / nrow(gwas)))
    cat(sprintf("  Output: %s\\n", out_file))
    cat(sprintf("\\nFirst 5 SNPs:\\n"))
    print(head(gwas_formatted, 5))
EOF
    
    if [ ! -f "${gwas_prefix}_formatted.txt" ]; then
        echo "❌ ERROR: Failed to format GWAS file" | tee -a format_gwas.log
        exit 1
    fi
    
    echo "" | tee -a format_gwas.log
    echo "✓ Formatting complete" | tee -a format_gwas.log
    
    mv format_gwas.log ${gwas_prefix}_format.log
    """
}
