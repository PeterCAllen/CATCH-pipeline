#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: catch_fdr.R <out.tsv> <catch_combined.tsv> [<catch_combined.tsv> ...]\n")
  cat("\n")
  cat("Applies Benjamini-Hochberg across all cell types AND traits within a dataset,\n")
  cat("which is the significance rule used in Li et al. A single pipeline run covers\n")
  cat("one trait, so this is run once over the outputs of all traits.\n")
  quit(status = 1)
}

out_file <- args[1]
in_files <- args[-1]

dt <- rbindlist(lapply(in_files, fread), use.names = TRUE, fill = TRUE)

for (col in c("group", "dataset", "trait", "CATCH_P")) {
  if (!col %in% names(dt)) {
    stop("Input is missing '", col, "'. Columns: ", paste(names(dt), collapse = ", "))
  }
}

if ("within_run_fdr" %in% names(dt)) dt[, within_run_fdr := NULL]

datasets <- unique(dt$dataset)
traits   <- unique(dt$trait)
cat("Datasets:", length(datasets), "| traits:", length(traits), "| rows:", nrow(dt), "\n")

if (length(traits) == 1) {
  cat("WARNING: only one trait present. The manuscript's FDR is computed across all\n")
  cat("traits within a dataset, so this will not reproduce the published thresholds.\n")
}

dt[, CATCH_FDR := p.adjust(CATCH_P, method = "fdr"), by = dataset]
dt[, significant := CATCH_FDR < 0.05]

setorder(dt, dataset, CATCH_FDR, CATCH_P)
fwrite(dt, out_file, sep = "\t")

cat("Wrote", out_file, "\n")
for (d in datasets) {
  sub <- dt[dataset == d]
  cat(sprintf("  %s: %d/%d trait-cell type pairs significant at 5%% FDR\n",
              d, sum(sub$significant, na.rm = TRUE), nrow(sub)))
}
