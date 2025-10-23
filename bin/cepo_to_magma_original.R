#!!!!!!!!!!! Be careful for MAGMA/top10.txt easy to be duplicated, only run once, if something, delete the txt output and then re-run

library(data.table)
library(dplyr)
library(tidyverse)
library(reshape2)

# Read gene coordinates file and process; get by Gene_coordinate_ready.py
gene_coordinates <- read_tsv("manual_analysis/geneMatrix.60048genes.tsv",
                             col_names = TRUE, col_types = 'cciicc') %>%
  mutate(start = ifelse(start - 100000 < 0, 0, start - 100000), end = end + 100000)


magma_top <- function(d, Cell_type, percent, total_genes) {
  n_genes_to_keep <- (total_genes * percent) %>% round()
  d_spe <- d %>% group_by(!!sym(Cell_type)) %>% top_n(., n_genes_to_keep, specificity)
  topNum <- 100 * percent
  d_spe %>% do(write_group_magma(., Cell_type, topNum))
}

write_group_magma <- function(df, Cell_type, topNum) {
  # Select only the cell type column and gene column
  df_subset <- select(df, all_of(Cell_type), gene)
  df_name <- make.names(unique(df_subset[[1]]))
  
  # Create output directory
  dir.create("manual_analysis/MAGMA/", showWarnings = FALSE)
  
  # Create output file path
  output_file <- paste0("manual_analysis/MAGMA/top", topNum, ".txt")
  
  # Prepare the data: transpose genes and add category name
  gene_list <- df_subset$gene
  result_df <- data.frame(Cat = df_name, t(gene_list), stringsAsFactors = FALSE)
  
  # Write to file (append mode)
  write_tsv(result_df, output_file, append = TRUE)
  return(df)
}


# Construct file path
file_path <- "results/cepo/cepo_stats.tsv"

# Check if file exists before processing

print(paste("Processing file:", file_path))

# Load data
dt <- fread(file_path)

# Use all gene types
gene_coordinates_filtered <- gene_coordinates %>%
    select(chr, start, end, Gene) %>%
    rename(gene = "Gene") %>%
    mutate(chr = paste0("chr", chr))

# Update column names to match gene coordinates
names(gene_coordinates_filtered) <- c("chr", "start", "end", "gene")
names(dt)[1]="gene"

# Merge data and gene coordinates
dt.exist <- inner_join(dt, gene_coordinates_filtered, by = "gene")

# remove gene_symbol
dt.exist <- dt.exist %>% select(-GENE_SYMBOL)

# Calculate number of genes and genes to keep
n_genes <- length(unique(dt.exist$gene))
n_genes_to_keep <- round(n_genes * 0.1)

# Melt data frame for processing
dt.melt <- melt(dt.exist, id.vars = c("gene", "chr", "start", "end"), 
                variable.name = "Cell_type", value.name = "specificity")
dt.melt <- dt.melt[, c(1, 5, 6, 2, 3, 4)]
dt.melt$Cell_type <- as.character(dt.melt$Cell_type)

# Define output directory and set it
output_dir <- paste0("manual_analysis/MAGMA/")

# Create the output directory if it does not exist
if (!dir.exists(output_dir)) {
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
print(paste("Created directory:", output_dir))
}

# Apply function to create MAGMA gene sets
dt.melt %>% magma_top("Cell_type", 0.1, n_genes)

print(paste("Completed processing for:", file_path))

