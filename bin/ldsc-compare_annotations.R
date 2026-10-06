#!/usr/bin/env Rscript

# Script to compare LDSC quantile results across cell types
# Usage: Rscript compare_annotations.R <results_dir> <gwas_name> <results_dir>

library(tidyverse)
library(ggrepel)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript compare_annotations.R <results_dir> <gwas_name>")
}

results_dir <- args[1]
gwas_name <- args[2]

cat("======================================\n")
cat("Comparing LDSC Quantile Results\n")
cat("======================================\n")
cat("Results directory:", results_dir, "\n")
cat("GWAS name:", gwas_name, "\n")
cat("Output directory:", results_dir, "\n\n")

# Function to extract results from a single file
extract_results <- function(file_path, annotation_name) {
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
  
  # Read results file
  results <- tryCatch({
    read.table(file_path, header = TRUE)
  }, error = function(e) {
    warning(paste("Error reading file:", file_path, "-", e$message))
    return(NULL)
  })
  
  if (is.null(results) || nrow(results) == 0) {
    warning(paste("No data in file:", file_path))
    return(NULL)
  }
  
  # Extract key metrics - overall enrichment and top quintile
  if (nrow(results) >= 5) {
    # Get top quintile (last row)
    top_quintile <- results[nrow(results), ]
    
    # Return summary data
    data.frame(
      annotation = annotation_name,
      # Overall metrics
      total_h2g = sum(results$h2g, na.rm = TRUE),
      # Top quintile metrics
      top_h2g = top_quintile$h2g,
      top_prop_h2g = top_quintile$prop_h2g,
      top_enr = top_quintile$enr,
      top_pval = top_quintile$enr_pval,
      stringsAsFactors = FALSE
    )
  } else {
    warning(paste("Not enough rows in results for", annotation_name))
    return(NULL)
  }
}

# Function to extract all quintile data from a file
extract_quintile_data <- function(file_path, annotation_name) {
  if (!file.exists(file_path)) {
    return(NULL)
  }
  
  results <- tryCatch({
    read.table(file_path, header = TRUE)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(results) || nrow(results) == 0) {
    return(NULL)
  }
  
  results$annotation <- annotation_name
  results$quintile <- 1:nrow(results)
  return(results)
}

# Find all .q5.txt files in the results directory
cat("Searching for quantile result files...\n")
q5_files <- list.files(results_dir, pattern = "\\.q5\\.txt$", full.names = TRUE)

if (length(q5_files) == 0) {
  stop("No .q5.txt files found in results directory: ", results_dir)
}

cat("Found", length(q5_files), "result files\n\n")

# Extract annotation names from filenames
annotation_data <- data.frame(
  file_path = q5_files,
  filename = basename(q5_files),
  stringsAsFactors = FALSE
) %>%
  mutate(
    annotation = sub("\\.q5\\.txt$", "", filename)
  )

cat("Processing annotations:\n")
print(annotation_data$annotation)
cat("\n")

# Process all annotations for summary statistics
all_results <- lapply(1:nrow(annotation_data), function(i) {
  extract_results(annotation_data$file_path[i], annotation_data$annotation[i])
}) %>%
  bind_rows()

if (nrow(all_results) == 0) {
  stop("No valid results extracted from files")
}

# Add significance indicators
all_results <- all_results %>%
  mutate(
    sig_level = case_when(
      top_pval < 1e-10 ~ "****",
      top_pval < 1e-8 ~ "***",
      top_pval < 1e-5 ~ "**",
      top_pval < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    # Log transform p-values for visualization
    neg_log_pval = -log10(top_pval)
  ) %>%
  # Sort by enrichment
  arrange(desc(top_enr))

# Write consolidated table to file
output_csv <- file.path(results_dir, paste0(gwas_name, "_annotation_comparison.csv"))
write.csv(all_results, file = output_csv, row.names = FALSE)
cat("Summary table written to:", output_csv, "\n")

# Generate comparison plots
cat("Generating comparison plots...\n")

# 1. Enrichment comparison
p1 <- ggplot(all_results, aes(x = reorder(annotation, top_enr), y = top_enr, fill = neg_log_pval)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sig_level), vjust = -0.5, size = 3) +
  scale_fill_viridis_c(name = "-log10(p-value)") +
  labs(title = paste(gwas_name, "- Annotation Enrichment Comparison"),
       x = "Cell Type", 
       y = "Enrichment (Top Quintile)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(file.path(results_dir, paste0(gwas_name, "_enrichment_comparison.pdf")), 
       p1, width = 12, height = 8)

# 2. Proportion of heritability explained
p2 <- ggplot(all_results, aes(x = reorder(annotation, top_prop_h2g), y = top_prop_h2g, fill = neg_log_pval)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_c(name = "-log10(p-value)") +
  labs(title = paste(gwas_name, "- Proportion of Heritability Explained"),
       x = "Cell Type", 
       y = "Proportion of h2g (Top Quintile)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(file.path(results_dir, paste0(gwas_name, "_prop_h2g_comparison.pdf")), 
       p2, width = 12, height = 8)

# 3. Scatter plot of enrichment vs proportion
p3 <- ggplot(all_results, aes(x = top_enr, y = top_prop_h2g, color = neg_log_pval)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = annotation), size = 3) +
  scale_color_viridis_c(name = "-log10(p-value)") +
  labs(title = paste(gwas_name, "- Enrichment vs Proportion of Heritability"),
       x = "Enrichment (Top Quintile)", 
       y = "Proportion of h2g (Top Quintile)") +
  theme_minimal()

ggsave(file.path(results_dir, paste0(gwas_name, "_enrichment_vs_prop.pdf")), 
       p3, width = 10, height = 8)

# Generate a heatmap for quintile patterns across annotations
cat("Generating quintile heatmap...\n")

quintile_data <- lapply(1:nrow(annotation_data), function(i) {
  extract_quintile_data(annotation_data$file_path[i], annotation_data$annotation[i])
}) %>%
  bind_rows()

if (nrow(quintile_data) > 0) {
  # Create heatmap of enrichment patterns
  p4 <- ggplot(quintile_data, aes(x = factor(quintile), y = annotation, fill = enr)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                        midpoint = 1, na.value = "gray90", name = "Enrichment") +
    labs(title = paste(gwas_name, "- Enrichment Patterns Across Quintiles"),
         x = "Quintile", 
         y = "Cell Type") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8))
  
  ggsave(file.path(results_dir, paste0(gwas_name, "_quintile_heatmap.pdf")), 
         p4, width = 12, height = max(10, nrow(all_results) * 0.3))
}

cat("\n======================================\n")
cat("Analysis complete!\n")
cat("Results saved to:", results_dir, "\n")
cat("======================================\n")