#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage: mbat_to_seismic_zstat.R <mbat_combined> <gene_coords> <out.tsv>\n")
  quit(status = 1)
}

mbat_file   <- args[1]
coords_file <- args[2]
out_file    <- args[3]

dt <- fread(mbat_file)

for (col in c("Gene", "P_mBATcombo")) {
  if (!col %in% names(dt)) {
    stop("ERROR: mBAT output missing '", col, "'. Columns: ", paste(names(dt), collapse = ", "))
  }
}

dt <- dt[!is.na(P_mBATcombo)]
if (nrow(dt) == 0) stop("ERROR: no genes with a non-missing P_mBATcombo.")

n_zero <- sum(dt$P_mBATcombo <= 0)
if (n_zero > 0) {
  pos <- dt$P_mBATcombo[dt$P_mBATcombo > 0]
  if (length(pos) == 0) stop("ERROR: every P_mBATcombo is zero.")
  floor_p <- min(pos) / 10
  cat("Flooring", n_zero, "zero p-values at", floor_p, "\n")
  dt[P_mBATcombo <= 0, P_mBATcombo := floor_p]
}

dt[P_mBATcombo > 1 - 1e-16, P_mBATcombo := 1 - 1e-16]

dt[, ZSTAT := qnorm(P_mBATcombo, lower.tail = FALSE)]

if (any(!is.finite(dt$ZSTAT))) {
  stop("ERROR: non-finite ZSTAT after transformation.")
}

coords <- fread(coords_file)
if (!all(c("Gene", "gene_name") %in% names(coords))) {
  stop("ERROR: gene coordinate file needs 'Gene' and 'gene_name'. Columns: ",
       paste(names(coords), collapse = ", "))
}

merged <- merge(dt[, .(Gene, ZSTAT)], coords[, .(Gene, gene_name)], by = "Gene", all.x = FALSE)

n_unmapped <- nrow(dt) - nrow(merged)
if (n_unmapped > 0) {
  cat("Dropped", n_unmapped, "genes with no symbol in the coordinate file\n")
}
if (nrow(merged) == 0) {
  stop("ERROR: no mBAT genes mapped to a gene symbol.")
}

merged <- merged[!is.na(gene_name) & nzchar(gene_name)]

out <- merged[, .(ZSTAT = max(ZSTAT)), by = .(GENE = gene_name)]

n_collapsed <- nrow(merged) - nrow(out)
if (n_collapsed > 0) {
  cat("Collapsed", n_collapsed, "duplicate gene symbols by max ZSTAT\n")
}

setorder(out, -ZSTAT)
fwrite(out, out_file, sep = "\t")

cat("Wrote", out_file, ":", nrow(out), "genes\n")
cat(sprintf("ZSTAT mean %.4f, sd %.4f, min %.3f, max %.3f\n",
            mean(out$ZSTAT), sd(out$ZSTAT), min(out$ZSTAT), max(out$ZSTAT)))
cat("Under the null ZSTAT should be ~N(0,1); large departures indicate a transform problem.\n")
