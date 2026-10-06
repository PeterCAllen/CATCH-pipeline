#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(argparse))

# --- Argument Parser ---
parser <- ArgumentParser(description = "Create MAGMA gene sets from Cepo results.")
parser$add_argument("--cepo_stats", type = "character", required = TRUE, 
                    help = "Path to the Cepo statistics TSV file (from run_cepo.R).")
parser$add_argument("--magma_gene_loc", type = "character", required = TRUE, 
                    help = "Path to MAGMA gene location file (Gene, chr, start, end, strand, gene_name) - NO HEADER.")
parser$add_argument("--top_pct", type = "double", default = 0.10, 
                    help = "Top percentage of genes to select for each cell type (default: 0.10 = 10%%).")
parser$add_argument("--output_file", type = "character", required = TRUE, 
                    help = "Path for the output MAGMA gene set file.")

args <- parser$parse_args()

# --- Read MAGMA gene location file (NO HEADER) ---
cat("Reading MAGMA gene location file...\n")
magma_genes <- read.table(args$magma_gene_loc, sep = "\t", header = FALSE, 
                          stringsAsFactors = FALSE, quote = "")
colnames(magma_genes) <- c("Gene", "chr", "start", "end", "strand", "gene_name")

cat(sprintf("Loaded %d genes from MAGMA gene location file\n", nrow(magma_genes)))

# --- Read Cepo stats ---
cat("Reading Cepo stats...\n")
cepo_sep <- if (grepl("\\.csv$", args$cepo_stats)) "," else "\t"
cepo_metrics <- read.table(args$cepo_stats, sep = cepo_sep, header = TRUE,
                           stringsAsFactors = FALSE, quote = "")

cat("Columns in Cepo stats:\n")
print(colnames(cepo_metrics))
cat(sprintf("Number of genes: %d\n", nrow(cepo_metrics)))

# Identify gene column
gene_col <- if ("gene" %in% colnames(cepo_metrics)) {
  "gene"
} else if ("Gene" %in% colnames(cepo_metrics)) {
  "Gene"
} else {
  stop("Could not find gene column in Cepo stats (expected 'gene' or 'Gene')")
}

# Get cell type columns (all columns except gene and gene_name columns)
# Exclude metadata columns: gene, Gene, gene_name, gene_id, gene_type, etc.
exclude_cols <- c("gene", "Gene", "gene_name", "gene_id", "gene_type", "chr", "start", "end")
cell_type_cols <- setdiff(colnames(cepo_metrics), exclude_cols)

cat(sprintf("Found %d cell types\n", length(cell_type_cols)))

# Pivot Cepo stats to long format
cepo_long <- cepo_metrics %>%
  pivot_longer(
    cols = all_of(cell_type_cols),
    names_to = "cell_type",
    values_to = "cepo_statistic"
  ) %>%
  rename(Gene = !!sym(gene_col))

# Filter genes that are in MAGMA gene location file
cat("Filtering genes present in MAGMA gene location file...\n")
cepo_long <- cepo_long %>%
  filter(Gene %in% magma_genes$Gene)

cat(sprintf("After filtering: %d gene-celltype pairs\n", nrow(cepo_long)))

# Select top genes per cell type
cat(sprintf("Selecting top %.1f%% of genes per cell type...\n", args$top_pct * 100))

top_genes <- cepo_long %>%
  group_by(cell_type) %>%
  slice_max(order_by = cepo_statistic, prop = args$top_pct, with_ties = FALSE) %>%
  ungroup()

cat(sprintf("Selected %d total gene-celltype pairs\n", nrow(top_genes)))

# Create gene sets (one per cell type)
magma_sets <- top_genes %>%
  group_by(cell_type) %>%
  summarise(genes = list(Gene), n_genes = n(), .groups = "drop")

cat(sprintf("\nGene set sizes per cell type:\n"))
print(magma_sets %>% select(cell_type, n_genes))

# Write to file in MAGMA format: cell_type\tgene1\tgene2\t...
cat(sprintf("\nWriting MAGMA gene set file to %s\n", args$output_file))

sink(args$output_file)
for (i in 1:nrow(magma_sets)) {
  cat(magma_sets$cell_type[i])
  cat("\t")
  cat(paste(magma_sets$genes[[i]], collapse = "\t"))
  cat("\n")
}
sink()

cat("MAGMA gene set creation complete.\n")
