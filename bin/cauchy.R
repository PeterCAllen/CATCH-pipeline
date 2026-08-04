#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

specs        <- character(0)
dataset_name <- NA_character_
gwas_prefix  <- NA_character_
out_prefix   <- NA_character_

i <- 1
while (i <= length(args)) {
  a <- args[i]
  if (a == "--method")       { specs <- c(specs, args[i + 1]); i <- i + 2 }
  else if (a == "--dataset") { dataset_name <- args[i + 1];    i <- i + 2 }
  else if (a == "--trait")   { gwas_prefix  <- args[i + 1];    i <- i + 2 }
  else if (a == "--out")     { out_prefix   <- args[i + 1];    i <- i + 2 }
  else stop("Unknown argument: ", a)
}

if (length(specs) == 0 || is.na(dataset_name) || is.na(gwas_prefix) || is.na(out_prefix)) {
  cat("Usage: cauchy.R --method NAME=FILE=KEYCOL=PCOL [--method ...] --dataset D --trait T --out PREFIX\n")
  quit(status = 1)
}

cauchy <- function(p) {
  if (all(is.na(p))) return(NA)
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
  if (cct.stat > 1e+15) (1 / cct.stat) / pi else 1 - pcauchy(cct.stat)
}

normalize_label <- function(x) {
  x <- gsub("[^A-Za-z0-9_-]", "_", as.character(x))
  x <- gsub("_+", "_", x)
  trimws(gsub("^[_-]+|[_-]+$", "", x))
}

parse_spec <- function(spec) {
  parts <- strsplit(spec, "=", fixed = TRUE)[[1]]
  if (length(parts) != 4) {
    stop("--method must be NAME=FILE=KEYCOL=PCOL, got: ", spec)
  }
  list(name = parts[1], file = parts[2], key = parts[3], pcol = parts[4])
}

methods <- lapply(specs, parse_spec)
method_names <- vapply(methods, function(m) m$name, character(1))

cat("Combining", length(methods), "methods:", paste(method_names, collapse = ", "), "\n\n")

read_method <- function(m) {
  if (!file.exists(m$file)) stop("File not found for method '", m$name, "': ", m$file)

  dt <- if (grepl("\\.csv$", m$file)) fread(m$file) else fread(m$file)

  if (!m$key %in% names(dt)) {
    stop("Method '", m$name, "': key column '", m$key, "' not in ", m$file,
         ". Columns: ", paste(names(dt), collapse = ", "))
  }
  if (!m$pcol %in% names(dt)) {
    stop("Method '", m$name, "': p-value column '", m$pcol, "' not in ", m$file,
         ". Columns: ", paste(names(dt), collapse = ", "))
  }

  out <- data.table(
    group = normalize_label(dt[[m$key]]),
    pval  = as.numeric(dt[[m$pcol]])
  )
  out <- out[!is.na(group) & nzchar(group)]

  n_dup <- sum(duplicated(out$group))
  if (n_dup > 0) {
    stop("Method '", m$name, "': ", n_dup, " duplicate cell types after label normalization.")
  }

  setnames(out, "pval", paste0(m$name, "_pval"))
  cat(sprintf("  %-24s %4d cell types from %s\n", m$name, nrow(out), basename(m$file)))
  out
}

tables <- lapply(methods, read_method)
cat("\n")

group_sets <- lapply(tables, function(t) t$group)
all_groups <- sort(Reduce(union, group_sets))
common     <- sort(Reduce(intersect, group_sets))

if (length(common) == 0) {
  cat("Cell types seen per method:\n")
  for (k in seq_along(methods)) {
    cat("  ", method_names[k], ": ", paste(utils::head(sort(group_sets[[k]]), 5), collapse = ", "), "\n", sep = "")
  }
  stop("No cell types are shared across all methods.")
}

if (length(common) < length(all_groups)) {
  cat("ERROR: cell type sets differ across methods.\n")
  for (k in seq_along(methods)) {
    missing <- setdiff(all_groups, group_sets[[k]])
    if (length(missing) > 0) {
      cat(sprintf("  %s is missing %d: %s\n", method_names[k], length(missing),
                  paste(missing, collapse = ", ")))
    }
  }
  stop("Refusing to silently drop cell types. Reconcile the labels upstream ",
       "(NORMALIZE_H5AD writes the canonical mapping) or investigate why a method dropped them.")
}

merged <- Reduce(function(x, y) merge(x, y, by = "group", all = FALSE), tables)
cat("Merged:", nrow(merged), "cell types across all", length(methods), "methods\n\n")

pval_cols <- paste0(method_names, "_pval")

for (col in pval_cols) {
  bad <- merged[[col]]
  if (any(is.na(bad))) cat("WARNING:", sum(is.na(bad)), "NA p-values in", col, "\n")
  if (any(!is.na(bad) & (bad < 0 | bad > 1))) stop("Out-of-range p-values in ", col)
}

merged[, CATCH_P := apply(.SD, 1, cauchy), .SDcols = pval_cols]
merged[, `:=`(dataset = dataset_name, trait = gwas_prefix)]

setcolorder(merged, c("group", "dataset", "trait", pval_cols, "CATCH_P"))
setorder(merged, CATCH_P)

merged[, within_run_fdr := p.adjust(CATCH_P, method = "fdr")]

out_tsv <- paste0(out_prefix, "_catch_combined.tsv")
fwrite(merged, out_tsv, sep = "\t")

plot_dt <- copy(merged)
plot_dt[, sig := cut(within_run_fdr, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                     labels = c("***", "**", "*", ""))]

p1 <- ggplot(plot_dt, aes(x = reorder(group, -CATCH_P), y = -log10(CATCH_P))) +
  geom_col(fill = "#3B6EA5", color = "black", linewidth = 0.2) +
  geom_text(aes(label = sig), hjust = -0.2, size = 3.5, na.rm = TRUE) +
  coord_flip() +
  theme_classic(base_size = 11) +
  labs(x = NULL,
       y = expression(-log[10]~"CATCH Cauchy P"),
       title = paste0("CATCH: ", gwas_prefix, " x ", dataset_name),
       caption = "Asterisks mark within-run FDR; the manuscript applies FDR across all traits (see catch_fdr.R)")

ggsave(paste0(out_prefix, "_catch_plot.png"), p1, dpi = 300,
       width = 7, height = max(3, 0.22 * nrow(plot_dt) + 1.5), limitsize = FALSE)

cat("Wrote", out_tsv, "\n")
cat("Cell types:", nrow(merged), "\n")
cat("Nominally significant (CATCH_P < 0.05):", sum(merged$CATCH_P < 0.05, na.rm = TRUE), "\n\n")
cat("Top 10 by CATCH_P:\n")
print(merged[1:min(10, nrow(merged)), c("group", pval_cols, "CATCH_P"), with = FALSE])
cat("\nNOTE: within_run_fdr covers this trait only. Li et al. control FDR at 5% across all\n")
cat("cell types AND traits within a dataset -- run bin/catch_fdr.R over all traits for that.\n")
