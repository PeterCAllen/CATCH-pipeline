// modules/mbat_to_tsv.nf
// Step 3: Convert mBAT results to TSV with Z-scores

process MBAT_TO_TSV {
    tag "${mbat_combined.simpleName}"
    label 'medium_mem'
    publishDir "${params.outdir}/scdrs/gwas", mode: 'copy'
    
    container "file://${projectDir}/environments/py-r-cepo-scdrs.sif"

    input:
    path mbat_combined

    output:
    path "*.tsv", emit: tsv_file
    path "*.log", emit: log

    script:
    def gwas_prefix = mbat_combined.simpleName.replaceAll(~/.gene.assoc.full$/, '')
    
    """
    echo "========================================" | tee mbat_to_tsv.log
    echo "Converting mBAT to TSV with Z-scores" | tee -a mbat_to_tsv.log
    echo "========================================" | tee -a mbat_to_tsv.log
    echo "Input: ${mbat_combined}" | tee -a mbat_to_tsv.log
    echo "Top genes: ${params.scdrs_top_genes}" | tee -a mbat_to_tsv.log
    echo "" | tee -a mbat_to_tsv.log
    
    Rscript - "${mbat_combined}" "${gwas_prefix}.tsv" "${params.scdrs_top_genes}" <<'RSCRIPT'
    suppressPackageStartupMessages({
      library(data.table)
      library(dplyr)
    })
    
    # ACATO function: Cauchy combination test for extreme p-values
    # More numerically stable than ACAT for very small p-values
    ACATO <- function(p) {
      if (all(is.na(p))) {
        return(NA)
      }
      p <- p[!is.na(p)]
      p[p == 1] <- 1 - 1e-16
      is.small <- (p < 1e-16)
      if (sum(is.small) == 0) {
        cct.stat <- sum(tan((0.5 - p) * pi)) / length(p)
      } else {
        cct.stat <- sum((1 / p[is.small]) / pi)
        cct.stat <- cct.stat + sum(tan((0.5 - p[!is.small]) * pi))
        cct.stat <- cct.stat / length(p)
      }
      if (cct.stat > 1e+15) {
        pval <- (1 / cct.stat) / pi
      } else {
        pval <- 1 - pcauchy(cct.stat)
      }
      pval
    }
    
    # Convert p-value to Z-score (two-tailed)
    pval_to_zscore <- function(p) {
      qnorm(p/2, lower.tail = FALSE) * sign(0.5 - p)
    }
    
    args <- commandArgs(trailingOnly = TRUE)
    mbat_file <- args[1]
    tsv_file <- args[2]
    top_n <- as.integer(args[3])
    
    dt <- fread(mbat_file)
    
    # Replace P = 0 with smallest non-zero to avoid numerical issues
    min_nonzero_p <- min(dt\$P_mBATcombo[dt\$P_mBATcombo > 0], na.rm = TRUE)
    dt[, P_mBATcombo := ifelse(P_mBATcombo == 0, min_nonzero_p/10, P_mBATcombo)]
    
    # Sort by P and take top N genes
    dt_sorted <- dt[order(P_mBATcombo)][1:min(.N, top_n)]
    
    # Convert each p-value to Z-score using ACATO for stability
    z_scores <- sapply(dt_sorted\$P_mBATcombo, function(p) {
      acato_p <- ACATO(p)
      pval_to_zscore(acato_p)
    })
    
    out <- data.table(GENE = dt_sorted\$Gene, Z = z_scores)
    
    fwrite(out, tsv_file, sep="\\t")
    
    cat(sprintf("✓ Created TSV with %d genes\\n", nrow(out)))
    cat(sprintf("  Min Z-score: %.3f\\n", min(out\$Z, na.rm=TRUE)))
    cat(sprintf("  Max Z-score: %.3f\\n", max(out\$Z, na.rm=TRUE)))
    cat(sprintf("  Mean Z-score: %.3f\\n", mean(out\$Z, na.rm=TRUE)))
    cat(sprintf("  Median Z-score: %.3f\\n", median(out\$Z, na.rm=TRUE)))
    
    # Check for extreme values
    n_inf <- sum(is.infinite(out\$Z))
    n_na <- sum(is.na(out\$Z))
    if (n_inf > 0) cat(sprintf("  ⚠️  Warning: %d Inf Z-scores\\n", n_inf))
    if (n_na > 0) cat(sprintf("  ⚠️  Warning: %d NA Z-scores\\n", n_na))
    
    # Report range of input p-values for diagnostics
    cat(sprintf("\\nInput p-value range:\\n"))
    cat(sprintf("  Min P-value: %.3e\\n", min(dt_sorted\$P_mBATcombo, na.rm=TRUE)))
    cat(sprintf("  Max P-value: %.3e\\n", max(dt_sorted\$P_mBATcombo, na.rm=TRUE)))
    RSCRIPT
    
    if [ ! -f "${gwas_prefix}.tsv" ]; then
        echo "❌ ERROR: Failed to create TSV file" | tee -a mbat_to_tsv.log
        exit 1
    fi
    
    echo "" | tee -a mbat_to_tsv.log
    echo "✓ TSV conversion complete" | tee -a mbat_to_tsv.log
    echo "  Output: ${gwas_prefix}.tsv" | tee -a mbat_to_tsv.log
    
    mv mbat_to_tsv.log ${gwas_prefix}_mbat_to_tsv.log
    """
}
