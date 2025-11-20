#!/usr/bin/env Rscript
# Combine P‑values from LDSC, MAGMA‑GSEA, scDRS via ACAT (Aggregated Cauchy Association Test)
# Updated to be generalizable for pipeline integration

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  cat("Usage: cauchy.R <ldsc_results> <magma_results> <scdrs_results> <h5ad_name> <gwas_prefix>\n")
  cat("  ldsc_results   : Path to LDSC annotation comparison CSV file\n")
  cat("  magma_results  : Path to MAGMA .gsa.out file\n")
  cat("  scdrs_results  : Path to scDRS .scdrs_group file\n")
  cat("  h5ad_name      : Name/identifier of the single-cell dataset\n")
  cat("  gwas_prefix    : GWAS trait identifier/prefix\n")
  quit(status = 1)
}

ldsc_file <- args[1]
magma_file <- args[2]
scdrs_file <- args[3]
dataset_name <- args[4]
gwas_prefix <- args[5]

cat("========================================\n")
cat("Cauchy Combination of P-values\n")
cat("========================================\n")
cat("LDSC results  :", ldsc_file, "\n")
cat("MAGMA results :", magma_file, "\n")
cat("scDRS results :", scdrs_file, "\n")
cat("Dataset       :", dataset_name, "\n")
cat("GWAS trait    :", gwas_prefix, "\n")
cat("\n")

# ---- cauchy Function ----
# Aggregated Cauchy Association Test for combining p-values
cauchy <- function(p) {
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

# ---- Read input files ----
cat("Reading input files...\n")

# LDSC results (annotation comparison CSV)
if (!file.exists(ldsc_file)) {
  stop("ERROR: LDSC results file not found: ", ldsc_file)
}
results_ldsc <- fread(ldsc_file)
cat("✓ Read LDSC results:", nrow(results_ldsc), "rows\n")

# MAGMA results (.gsa.out)
if (!file.exists(magma_file)) {
  stop("ERROR: MAGMA results file not found: ", magma_file)
}
results_magma <- fread(magma_file, skip = 3)
cat("✓ Read MAGMA results:", nrow(results_magma), "rows\n")

# scDRS results (.scdrs_group)
if (!file.exists(scdrs_file)) {
  stop("ERROR: scDRS results file not found: ", scdrs_file)
}
results_scdrs <- fread(scdrs_file)
cat("✓ Read scDRS results:", nrow(results_scdrs), "rows\n")
cat("\n")

# ---- Standardize column names ----
cat("Standardizing column names...\n")

# LDSC: rename "annotation" to "group"
if ("annotation" %in% names(results_ldsc)) {
  setnames(results_ldsc, "annotation", "group")
}

# MAGMA: rename "VARIABLE" to "group"
if ("VARIABLE" %in% names(results_magma)) {
  setnames(results_magma, "VARIABLE", "group")
}

# scDRS: should already have "group" column
if (!"group" %in% names(results_scdrs)) {
  stop("ERROR: scDRS results missing 'group' column")
}

# ---- Extract relevant p-value columns ----
cat("Extracting p-values...\n")

# LDSC: use "top_pval" (the best quintile p-value)
if (!"top_pval" %in% names(results_ldsc)) {
  stop("ERROR: LDSC results missing 'top_pval' column")
}
ldsc_cols <- results_ldsc[, .(group, ldsc_pval = top_pval)]

# MAGMA: use "P" (gene-set association p-value)
if (!"P" %in% names(results_magma)) {
  stop("ERROR: MAGMA results missing 'P' column")
}
magma_cols <- results_magma[, .(group, magma_pval = P)]

# scDRS: use "assoc_mcp" (MC-corrected cell type association p-value)
if (!"assoc_mcp" %in% names(results_scdrs)) {
  stop("ERROR: scDRS results missing 'assoc_mcp' column")
}
scdrs_cols <- results_scdrs[, .(group, scdrs_pval = assoc_mcp)]

cat("✓ LDSC p-values extracted for", nrow(ldsc_cols), "cell types\n")
cat("✓ MAGMA p-values extracted for", nrow(magma_cols), "cell types\n")
cat("✓ scDRS p-values extracted for", nrow(scdrs_cols), "cell types\n")
cat("\n")

# ---- Merge all three data.tables on cell type/group ----
cat("Merging results by cell type...\n")

# Use inner join (all = FALSE) to keep only cell types present in all three analyses
merged_data <- merge(ldsc_cols, magma_cols, by = "group", all = FALSE)
merged_data <- merge(merged_data, scdrs_cols, by = "group", all = FALSE)

cat("✓ Merged data:", nrow(merged_data), "cell types present in all three analyses\n")

if (nrow(merged_data) == 0) {
  stop("ERROR: No common cell types found across all three analyses")
}
cat("\n")

# ---- Calculate cauchy combined p-value ----
cat("Computing Cauchy combined p-values...\n")

merged_data[, CAUCHY_P := {
  pvals <- c(ldsc_pval, magma_pval, scdrs_pval)
  cauchy(pvals)
}, by = 1:nrow(merged_data)]

cat("✓ Computed cauchy p-values\n")
cat("\n")

# ---- Add metadata columns ----
merged_data[, `:=`(dataset = dataset_name, trait = gwas_prefix)]

# Reorder columns for clarity
setcolorder(merged_data, c("group", "dataset", "trait", "ldsc_pval", "magma_pval", "scdrs_pval", "CAUCHY_P"))


# ---- Plot cauchy p-values ----
merged_data <- merged_data %>%
    mutate(assoc_fdr_0.05 = p.adjust(CAUCHY_P, method = "fdr"),
           assoc_fdr_0.05_fig = cut(
            assoc_fdr_0.05,
            breaks = c(-Inf, 0.001, 0.05, 0.1, 1),
            labels = c("***", "**", "*", "")
    )
)

p1 <- ggplot(merged_data, aes(x = group, y = -log10(assoc_fdr_0.05), fill = group)) +
  geom_col(color = "black") +
  geom_text(aes(label = assoc_fdr_0.05_fig),
            vjust = -0.3, size = 4, na.rm = TRUE) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(face = "bold", angle = 90, hjust = 1),
    panel.border = element_rect(fill = NA, color = "black", size = 1),
    legend.position = "none"
  ) +
  labs(
    y = expression(-log[10]~"FDR-adjusted ACAT-O P"),
    title = "Significant associations across cell types"
  )

ggsave(
  filename = paste0(gwas_prefix, "_cauchy_cauchy_fdr_plot.png"),
  plot = p1,
  dpi = 300
)

# ---- Write output ----
output_file <- paste0(gwas_prefix, "_cauchy_combined.tsv")
fwrite(merged_data, file = output_file, sep = "\t")

cat("========================================\n")
cat("Results Summary\n")
cat("========================================\n")
cat("Cell types analyzed:", nrow(merged_data), "\n")
cat("Significant (CAUCHY_P < 0.05):", sum(merged_data$CAUCHY_P < 0.05, na.rm = TRUE), "\n")
cat("\n")
cat("Top 10 cell types by cauchy p-value:\n")
top_results <- merged_data[order(CAUCHY_P)][1:min(10, nrow(merged_data))]
print(top_results[, .(group, ldsc_pval, magma_pval, scdrs_pval, CAUCHY_P)])
cat("\n")
cat("✓ Saved combined results to:", output_file, "\n")
cat("========================================\n")