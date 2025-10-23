// modules/combine_mbat.nf
// Step 2: Combine mBAT results from all chromosomes

process COMBINE_MBAT {
    tag "${gwas_prefix}"
    label 'medium_mem'
    publishDir "${params.outdir}/scdrs/gwas", mode: 'copy'
    
    container "file://${projectDir}/environments/py-r-cepo-scdrs.sif"

    input:
    tuple val(gwas_prefix), path(mbat_files)

    output:
    path "*.gene.assoc.full.mbat", emit: mbat_combined
    path "*.log", emit: log

    script:
    """
    echo "========================================" | tee combine_mbat.log
    echo "Combining mBAT results from 22 chromosomes" | tee -a combine_mbat.log
    echo "========================================" | tee -a combine_mbat.log
    echo "GWAS: ${gwas_prefix}" | tee -a combine_mbat.log
    echo "Input files: ${mbat_files.size()}" | tee -a combine_mbat.log
    echo "" | tee -a combine_mbat.log
    
    Rscript - "${gwas_prefix}" <<'RSCRIPT'
    library(tidyverse)
    
    args <- commandArgs(trailingOnly = TRUE)
    gwas_prefix <- args[1]
    
    # List all mBAT files
    file_paths <- list.files(pattern = "\\\\.gene\\\\.assoc\\\\.mbat\$", full.names = TRUE)
    
    cat(sprintf("Found %d mBAT files\\n", length(file_paths)))
    
    # Read and combine all chromosome files
    gwas_combined <- map_dfr(file_paths, read_delim, delim = "\\t", show_col_types = FALSE)
    
    out_file <- paste0(gwas_prefix, ".gene.assoc.full.mbat")
    write_delim(gwas_combined, out_file, delim = "\\t")
    
    cat(sprintf("✓ Combined %d genes from %d chromosomes\\n", nrow(gwas_combined), length(file_paths)))
    cat(sprintf("  Output: %s\\n", out_file))
    cat(sprintf("\\nFirst 5 genes:\\n"))
    print(head(gwas_combined, 5))
    RSCRIPT
    
    if [ ! -f "${gwas_prefix}.gene.assoc.full.mbat" ]; then
        echo "❌ ERROR: Failed to combine mBAT results" | tee -a combine_mbat.log
        exit 1
    fi
    
    echo "" | tee -a combine_mbat.log
    echo "✓ Combining complete" | tee -a combine_mbat.log
    
    mv combine_mbat.log ${gwas_prefix}_combine_mbat.log
    """
}
