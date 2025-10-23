#!/usr/bin/env Rscript

# Cepo analysis script for conda environment
# Performs cell type prioritization analysis on single-cell data

# Load additional required libraries

suppressPackageStartupMessages({
  library(Cepo)
  library(anndataR)
  library(SingleCellExperiment)
})

print("=== Cepo Analysis ===")

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Parse arguments manually
input_h5ad <- NULL
cell_type_col <- NULL
output_rds <- NULL
output_tsv <- NULL
gene_coords <- NULL

for (i in seq_along(args)) {
  if (args[i] == "--input_h5ad") input_h5ad <- args[i+1]
  if (args[i] == "--cell_type_col") cell_type_col <- args[i+1]
  if (args[i] == "--output_rds") output_rds <- args[i+1]
  if (args[i] == "--output_tsv") output_tsv <- args[i+1]
  if (args[i] == "--gene_coords") gene_coords <- args[i+1]
}

print(paste("Input file:", input_h5ad))
print(paste("Cell type column:", cell_type_col))

# Read h5ad file
print("Reading h5ad file with anndataR...")
sce <- read_h5ad(input_h5ad, as = "SingleCellExperiment")

print(paste("Data loaded:", ncol(sce), "cells x", nrow(sce), "genes"))

# Check if cell type column exists
if (!cell_type_col %in% colnames(colData(sce))) {
  available_cols <- colnames(colData(sce))
  print(paste("Available columns:", paste(available_cols, collapse = ", ")))
  stop(paste("Cell type column", cell_type_col, "not found"))
}

cell_types <- colData(sce)[[cell_type_col]]
print(paste("Cell types found:", length(unique(cell_types))))
print("Cell type distribution:")
print(table(cell_types))

print("Running Cepo analysis...")

# Convert to more compatible format for Cepo
print("Converting matrix format for compatibility...")
print("Available assays:")
print(assayNames(sce))
print("Class of main assay:")
print(class(assay(sce, "X")))

# Convert to a dense matrix that Cepo can handle
print("Converting to dense matrix...")
counts_matrix <- as.matrix(assay(sce, "X"))
print(paste("Matrix dimensions:", nrow(counts_matrix), "x", ncol(counts_matrix)))
print(paste("Matrix class:", class(counts_matrix)))

# Try passing just the matrix directly to Cepo instead of SCE
print("Running Cepo with matrix directly...")
# Run Cepo analysis with just the matrix
cepo_results <- Cepo(
  exprsMat = counts_matrix,
  cellTypes = cell_types,
  prefilter_scExp = FALSE
)

print("Saving results...")
# Save results
saveRDS(cepo_results, file = output_rds)

# Extract and save summary statistics using base R
cepo_stats <- cepo_results$stats
cepo_stats$gene <- rownames(cepo_stats)

# Add gene symbol mapping if gene coordinates file is provided
if (!is.null(gene_coords) && file.exists(gene_coords)) {
  print("Loading gene coordinate mapping...")
  gene_map <- read.table(gene_coords, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  # Merge with cepo_stats to add gene symbols
  cepo_stats <- merge(cepo_stats, gene_map[, c("Gene", "gene_name")], 
                      by.x = "gene", by.y = "Gene", all.x = TRUE)

  # Reorder columns to put gene and gene_name first
  cepo_stats <- cepo_stats[, c("gene", "gene_name", setdiff(names(cepo_stats), c("gene", "gene_name")))]
  
  print("Added gene symbol mapping")
} else {
  # Reorder columns using base R (original behavior)
  cepo_stats <- cepo_stats[, c("gene", setdiff(names(cepo_stats), "gene"))]
}

# Write using base R
write.table(cepo_stats, output_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

print("=== Cepo analysis completed successfully ===")
print(paste("Results saved to:", output_rds))
print(paste("Statistics saved to:", output_tsv))
