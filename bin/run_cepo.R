#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Cepo)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
})

args <- commandArgs(trailingOnly = TRUE)

input_h5ad      <- NULL
cell_type_col   <- NULL
output_rds      <- NULL
output_stats    <- NULL
output_pvalues  <- NULL
compute_pvalue  <- 100
prefilter_pzero <- 0.4

for (i in seq_along(args)) {
  if (args[i] == "--input_h5ad")      input_h5ad      <- args[i + 1]
  if (args[i] == "--cell_type_col")   cell_type_col   <- args[i + 1]
  if (args[i] == "--output_rds")      output_rds      <- args[i + 1]
  if (args[i] == "--output_stats")    output_stats    <- args[i + 1]
  if (args[i] == "--output_pvalues")  output_pvalues  <- args[i + 1]
  if (args[i] == "--compute_pvalue")  compute_pvalue  <- as.integer(args[i + 1])
  if (args[i] == "--prefilter_pzero") prefilter_pzero <- as.numeric(args[i + 1])
}

cat("Input h5ad      :", input_h5ad, "\n")
cat("Cell type column:", cell_type_col, "\n")
cat("computePvalue   :", compute_pvalue, "\n")
cat("prefilter_pzero :", prefilter_pzero, "\n")

sce <- readH5AD(input_h5ad, reader = "R")
cat("Loaded", ncol(sce), "cells x", nrow(sce), "genes\n")

if (!cell_type_col %in% colnames(colData(sce))) {
  stop("Cell type column '", cell_type_col, "' not found. Available: ",
       paste(colnames(colData(sce)), collapse = ", "))
}

assay_name <- if ("X" %in% assayNames(sce)) "X" else assayNames(sce)[1]
cat("Using assay:", assay_name, "\n")

exprs_mat <- as.matrix(assay(sce, assay_name))
cell_types <- as.character(colData(sce)[[cell_type_col]])

cat("Cell type distribution:\n")
print(table(cell_types))

set.seed(602)

cepo_results <- Cepo(
  exprsMat        = exprs_mat,
  cellTypes       = cell_types,
  computePvalue   = compute_pvalue,
  prefilter_pzero = prefilter_pzero
)

saveRDS(cepo_results, file = output_rds)

write_matrix <- function(mat, path) {
  df <- as.data.frame(mat)
  df$gene <- rownames(df)
  df <- df[, c("gene", setdiff(names(df), "gene")), drop = FALSE]
  rownames(df) <- NULL
  write.table(df, path, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
  cat("Wrote", path, ":", nrow(df), "genes x", ncol(df) - 1, "cell types\n")
}

write_matrix(cepo_results$stats, output_stats)

if (is.null(cepo_results$pvalues)) {
  stop("Cepo returned no p-value matrix; computePvalue must be > 0 for the cepo_s variant.")
}
write_matrix(cepo_results$pvalues, output_pvalues)

cat("Cepo complete.\n")
