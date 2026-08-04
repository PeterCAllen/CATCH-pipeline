#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: munge_gwas_for_ldsc.R <in_with_rsids.txt> <out_for_ldsc.txt> [N]\n")
  quit(status = 1)
}

in_file  <- args[1]
out_file <- args[2]
n_arg    <- if (length(args) >= 3 && nzchar(args[3])) as.numeric(args[3]) else NA_real_

ALIASES <- list(
  SNP  = c("SNP", "RSID", "RS", "SNPID", "MARKERNAME", "VARIANT_ID", "ID"),
  A1   = c("A1", "EFFECT_ALLELE", "ALLELE1", "EA", "TESTED_ALLELE"),
  A2   = c("A2", "OTHER_ALLELE", "ALLELE2", "NEA", "REF", "NON_EFFECT_ALLELE"),
  N    = c("N", "NEFF", "N_TOTAL", "SAMPLESIZE", "TOTALSAMPLESIZE"),
  Z    = c("Z", "ZSCORE", "Z_SCORE", "ZSTAT"),
  BETA = c("BETA", "B", "EFFECT", "EFFECTS", "LOG_ODDS"),
  OR   = c("OR", "ODDS_RATIO"),
  SE   = c("SE", "STDERR", "STANDARD_ERROR"),
  P    = c("P", "PVAL", "P_VALUE", "PVALUE", "P_BOLT_LMM", "P_LINREG"),
  FRQ  = c("FRQ", "FREQ", "MAF", "EAF", "A1FREQ", "EFFECT_ALLELE_FREQUENCY"),
  INFO = c("INFO", "IMPINFO", "RSQ")
)

dt <- fread(in_file, header = TRUE, data.table = TRUE)
orig <- toupper(gsub("[^A-Za-z0-9_]", "_", names(dt)))

resolved <- character(0)
for (canon in names(ALIASES)) {
  hit <- which(orig %in% ALIASES[[canon]])
  if (length(hit) > 0) {
    setnames(dt, names(dt)[hit[1]], canon)
    resolved <- c(resolved, canon)
  }
}

cat("Resolved columns:", paste(resolved, collapse = ", "), "\n")

for (req in c("SNP", "A1", "A2")) {
  if (!req %in% names(dt)) {
    stop("ERROR: required column '", req, "' could not be resolved from header: ",
         paste(names(dt), collapse = ", "))
  }
}

dt[, A1 := toupper(as.character(A1))]
dt[, A2 := toupper(as.character(A2))]

if (!"Z" %in% names(dt)) {
  if (all(c("BETA", "SE") %in% names(dt))) {
    cat("Deriving Z from BETA/SE\n")
    dt[, Z := BETA / SE]
  } else if (all(c("OR", "SE") %in% names(dt))) {
    cat("Deriving Z from log(OR)/SE\n")
    dt[, Z := log(OR) / SE]
  }
}

if (!"P" %in% names(dt)) {
  if ("Z" %in% names(dt)) {
    cat("Deriving P from Z\n")
    dt[, P := 2 * pnorm(-abs(Z))]
  } else {
    stop("ERROR: no P column and no way to derive one (need Z, or BETA/SE, or OR/SE).")
  }
}

if (!"N" %in% names(dt)) {
  if (is.na(n_arg)) {
    stop("ERROR: no N column in the GWAS and no N supplied on the command line.")
  }
  cat("Adding constant N =", n_arg, "\n")
  dt[, N := n_arg]
}

if (!"Z" %in% names(dt) && !any(c("BETA", "OR") %in% names(dt))) {
  stop("ERROR: no signed effect column (Z, BETA or OR); munge_sumstats.py cannot run.")
}

keep <- intersect(c("SNP", "A1", "A2", "N", "Z", "BETA", "OR", "SE", "P", "FRQ", "INFO"), names(dt))
out <- dt[, ..keep]

out <- out[!is.na(SNP) & nzchar(SNP) & SNP != "."]
out <- out[!duplicated(SNP)]
out <- out[!is.na(P) & P > 0 & P <= 1]

fwrite(out, out_file, sep = "\t", quote = FALSE, na = "NA")

cat("Wrote", nrow(out), "variants to", out_file, "\n")
cat("Columns:", paste(names(out), collapse = ", "), "\n")
