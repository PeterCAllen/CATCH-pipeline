#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage: esw_star.R <stats.csv> <pvalues.csv> <out.csv>\n")
  quit(status = 1)
}

stats_file <- args[1]
pval_file  <- args[2]
out_file   <- args[3]

esw_star <- function(esw, pvals) {
  pval_mask    <- (pvals <= 0.05)
  binzero_mask <- (esw > 0)
  mask <- !(pval_mask & binzero_mask)

  esw_nominal <- esw
  esw_nominal[mask] <- NA

  esw_ranked <- apply(esw_nominal, 2, rank, na.last = "keep", ties.method = "average")
  esw_ranked[is.na(esw_ranked)] <- 0

  col_max <- apply(esw_ranked, 2, max)
  col_max[col_max == 0] <- 1

  as.data.frame(sweep(esw_ranked, 2, col_max, FUN = "/"))
}

dt_stats <- fread(stats_file)
dt_pvals <- fread(pval_file)

setnames(dt_stats, 1, "gene")
setnames(dt_pvals, 1, "gene")

common_genes <- intersect(dt_stats$gene, dt_pvals$gene)
if (length(common_genes) == 0) {
  stop("No genes shared between the Cepo stats and p-value matrices.")
}
if (length(common_genes) < nrow(dt_stats)) {
  cat("Restricting to", length(common_genes), "genes present in both matrices\n")
}

dt_stats <- dt_stats[match(common_genes, gene)]
dt_pvals <- dt_pvals[match(common_genes, gene)]

annot_cols <- intersect(setdiff(names(dt_stats), "gene"), setdiff(names(dt_pvals), "gene"))
if (length(annot_cols) == 0) {
  stop("No shared cell type columns between the Cepo stats and p-value matrices.")
}

esw   <- as.matrix(dt_stats[, ..annot_cols])
pvals <- as.matrix(dt_pvals[, ..annot_cols])
rownames(esw) <- common_genes

out <- esw_star(esw, pvals)
out <- cbind(data.frame(gene = common_genes, stringsAsFactors = FALSE), out)

fwrite(out, out_file, sep = ",", quote = FALSE)

cat("Wrote", out_file, ":", nrow(out), "genes x", length(annot_cols), "cell types\n")
