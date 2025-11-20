#!/usr/bin/env python
"""
Generate genome-build-specific gene coordinate files from geneMatrix.tsv.gz
Based on manual analysis scripts: 0-Gene_coordinate_ready-hg19.py and 0-Gene_coordinate_ready-hg38.py
"""
import pandas as pd
import argparse
import time

def log_message(message):
    """Prints a formatted log message with a timestamp."""
    timestamp = time.strftime('%Y-%m-%d %H:%M:%S')
    print(f"--- {timestamp} --- {message}", flush=True)

def main(args):
    """
    Parses geneMatrix.tsv.gz to create hg19 coordinate files for MAGMA, LDSC, CEPO, and scDRS.
    
    NOTE: ALL analyses (MAGMA, LDSC, scDRS/mBAT) use hg19 coordinates for reference matching.
    The genome_build parameter is used for:
    1. Logging purposes
    2. SNP matching in GWAS files (SNPs use their original coordinates)
    3. Gene coordinates are ALWAYS hg19 (for reference panel matching)
    """
    log_message("Python script started.")
    log_message(f"Processing for genome build: {args.genome_build}")
    
    # Load gene matrix file 
    log_message(f"Reading gene matrix file: {args.gene_matrix}")
    read_start_time = time.time()
    
    gene_coordinates = pd.read_csv(args.gene_matrix, delimiter='\t', compression='gzip')
    
    read_end_time = time.time()
    log_message(f"Finished reading gene matrix in {read_end_time - read_start_time:.2f} seconds.")
    log_message(f"Initial dataframe shape: {gene_coordinates.shape}")

    # Validate genome build parameter
    if args.genome_build not in ['hg19', 'GRCh37', 'hg38', 'GRCh38']:
        raise ValueError(f"Unknown genome build: {args.genome_build}. Use hg19, GRCh37, hg38, or GRCh38.")
    
    # ALL analyses use hg19 coordinates (for reference file matching)
    # Note: GWAS SNPs retain their original coordinates (hg19 or hg38)
    coord_chr_col = 'hg19g0'
    coord_start_col = 'g1'
    coord_end_col = 'g2'
    coord_strand_col = 'gstr'
    
    log_message(f"All analyses will use hg19 gene coordinates (regardless of GWAS build)")
    log_message(f"GWAS SNPs retain their original {args.genome_build} coordinates")
    
    # Sort the DataFrame and drop duplicates
    gene_coordinates.drop_duplicates(inplace=True)
    log_message(f"After dropping duplicates: {gene_coordinates.shape}")
    
    # Remove rows with NaN in hg19 coordinates
    gene_coordinates.dropna(subset=[coord_start_col, coord_end_col], inplace=True)
    log_message(f"After removing NaN coordinates: {gene_coordinates.shape}")
    
    # Prepare gene coordinates (always hg19 for all analyses)
    df_coords = gene_coordinates.copy()
    df_coords.rename(columns={coord_chr_col: 'chr', 'ensgid': 'Gene'}, inplace=True)
    df_coords['start'] = df_coords[coord_start_col].astype(int)
    df_coords['end'] = df_coords[coord_end_col].astype(int)
    df_coords['chr'] = df_coords['chr'].str.replace('chr', '', regex=False)
    df_coords['chr'] = pd.Categorical(
        df_coords['chr'], 
        categories=[str(i) for i in range(1, 23)] + ['X', 'Y', 'M'],
        ordered=True
    )
    df_coords = df_coords[~df_coords['chr'].isin(['M', 'X', 'Y'])]
    df_coords.sort_values(by=['chr', 'start', 'end'], inplace=True)
    log_message(f"Gene coordinates (hg19): {df_coords.shape}")
    
    # === Output 1: LDSC/mBAT version (chr, start, end, Gene) - NO HEADER - hg19 ===
    log_message("Creating LDSC/mBAT gene coordinate file (hg19)...")
    df_ldsc_out = df_coords[['chr', 'start', 'end', 'Gene']].copy()
    df_ldsc_out.to_csv(args.output_ldsc, sep='\t', index=False, header=False)
    log_message(f"LDSC gene coordinate file saved to {args.output_ldsc} ({len(df_ldsc_out)} genes, hg19 coordinates)")
    
    # === Output 2: MAGMA version (Gene, chr, start, end, strand, gene_name) - NO HEADER - hg19 ===
    log_message("Creating MAGMA gene coordinate file (hg19)...")
    df_magma_out = df_coords[['Gene', 'chr', 'start', 'end', coord_strand_col, 'gene_name']].copy()
    df_magma_out.to_csv(args.output_magma, sep='\t', index=False, header=False)
    log_message(f"MAGMA gene coordinate file saved to {args.output_magma} ({len(df_magma_out)} genes, hg19 coordinates)")
    
    # === Output 3: CEPO/scDRS version (Gene, chr, start, end, gene_type, gene_name) - WITH HEADER - hg19 ===
    # Note: CEPO analysis works with gene expression and doesn't need specific genomic coordinates
    # scDRS uses mBAT which requires hg19 reference matching, so we use hg19 coordinates
    log_message("Creating CEPO/scDRS gene coordinate file (hg19)...")
    df_cepo_out = df_coords[['Gene', 'chr', 'start', 'end', 'gene_type', 'gene_name']].copy()
    df_cepo_out.to_csv(args.output_cepo, sep='\t', index=False, header=True)
    log_message(f"CEPO gene coordinate file saved to {args.output_cepo} ({len(df_cepo_out)} genes, hg19 coordinates)")
    
    log_message("Python script finished successfully.")
    log_message("Summary:")
    log_message("  - ALL gene coordinates use hg19 (for reference file matching)")
    log_message(f"  - GWAS build ({args.genome_build}) is preserved for SNP coordinates")
    log_message("  - This ensures all analyses use compatible hg19 reference panels")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Generate genome-build-specific gene coordinate files from geneMatrix.tsv.gz"
    )
    parser.add_argument('--gene_matrix', type=str, required=True, 
                        help='Path to the input gene matrix file (geneMatrix.tsv.gz)')
    parser.add_argument('--genome_build', type=str, required=True, 
                        choices=['hg19', 'GRCh37', 'hg38', 'GRCh38'],
                        help='Genome build (hg19, GRCh37, hg38, or GRCh38)')
    parser.add_argument('--output_magma', type=str, required=True, 
                        help='Path for the output MAGMA gene location file (no header)')
    parser.add_argument('--output_ldsc', type=str, required=True, 
                        help='Path for the output LDSC gene location file (no header)')
    parser.add_argument('--output_cepo', type=str, required=True, 
                        help='Path for the output CEPO gene coordinate file (with header)')
    
    args = parser.parse_args()
    main(args)