#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(seismicGWAS)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

input_h5ad    <- NULL
cell_type_col <- NULL
zstat_file    <- NULL
output_tsv    <- NULL
output_sscore <- NULL
assay_name    <- "logcounts"
lognorm       <- FALSE

for (i in seq_along(args)) {
  if (args[i] == "--input_h5ad")    input_h5ad    <- args[i + 1]
  if (args[i] == "--cell_type_col") cell_type_col <- args[i + 1]
  if (args[i] == "--zstat_file")    zstat_file    <- args[i + 1]
  if (args[i] == "--output_tsv")    output_tsv    <- args[i + 1]
  if (args[i] == "--output_sscore") output_sscore <- args[i + 1]
  if (args[i] == "--assay_name")    assay_name    <- args[i + 1]
  if (args[i] == "--lognorm")       lognorm       <- tolower(args[i + 1]) %in% c("true", "yes", "1")
}

sce <- readH5AD(input_h5ad, reader = "R")
cat("Loaded", ncol(sce), "cells x", nrow(sce), "genes\n")

if (!cell_type_col %in% colnames(colData(sce))) {
  stop("Cell type column '", cell_type_col, "' not found. Available: ",
       paste(colnames(colData(sce)), collapse = ", "))
}

if (lognorm) {
  cat("Running scater::logNormCounts\n")
  if (!"counts" %in% assayNames(sce)) {
    assay(sce, "counts") <- assay(sce, assayNames(sce)[1])
  }
  sce <- scater::logNormCounts(sce)
  assay_name <- "logcounts"
} else if (!assay_name %in% assayNames(sce)) {
  src <- if ("X" %in% assayNames(sce)) "X" else assayNames(sce)[1]
  cat("Assay '", assay_name, "' absent; aliasing '", src, "' to it\n", sep = "")
  assay(sce, assay_name) <- assay(sce, src)
}

cat("Computing seismic specificity on assay:", assay_name, "\n")
sscore <- calc_specificity(sce, assay_name = assay_name, ct_label_col = cell_type_col)
cat("Specificity matrix:", nrow(sscore), "genes x", ncol(sscore), "cell types\n")

if (!is.null(output_sscore)) {
  ss <- as.data.frame(as.matrix(sscore))
  ss$gene <- rownames(ss)
  ss <- ss[, c("gene", setdiff(names(ss), "gene")), drop = FALSE]
  fwrite(ss, output_sscore, sep = ",")
}

mbat <- fread(zstat_file)
if (!all(c("GENE", "ZSTAT") %in% names(mbat))) {
  stop("ERROR: z-statistic file needs 'GENE' and 'ZSTAT'. Columns: ",
       paste(names(mbat), collapse = ", "))
}

overlap <- length(intersect(rownames(sscore), mbat$GENE))
frac <- overlap / nrow(sscore)
cat(sprintf("Gene overlap: %d/%d specificity genes present in mBAT (%.1f%%)\n",
            overlap, nrow(sscore), 100 * frac))
if (overlap == 0) {
  stop("ERROR: no genes shared between the specificity matrix and the mBAT z-statistics. ",
       "Check that both sides use the same gene identifier type.")
}
if (frac < 0.8) {
  cat("WARNING: overlap below seismic's 80% threshold; results may be unreliable.\n")
}

res <- get_ct_trait_associations(
  sscore,
  mbat,
  magma_gene_col = "GENE",
  magma_z_col    = "ZSTAT"
)

res <- as.data.table(res)
fwrite(res, output_tsv, sep = "\t")

cat("Wrote", output_tsv, ":", nrow(res), "cell types\n")
print(utils::head(res[order(pvalue)], 10))
