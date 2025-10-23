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
    Parses geneMatrix.tsv.gz to create genome-build-specific coordinate files for MAGMA, LDSC, and CEPO.
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

    # Select columns based on genome build
    if args.genome_build in ['hg19', 'GRCh37']:
        chr_col = 'hg19g0'
        start_col = 'g1'
        end_col = 'g2'
        strand_col = 'gstr'
    elif args.genome_build in ['hg38', 'GRCh38']:
        chr_col = 'hg38h0'  
        start_col = 'h1'
        end_col = 'h2'
        strand_col = 'hstr'
    else:
        raise ValueError(f"Unknown genome build: {args.genome_build}. Use hg19, GRCh37, hg38, or GRCh38.")
    
    # Sort the DataFrame based on start position and drop duplicates
    gene_coordinates.sort_values(by=start_col, inplace=True)
    gene_coordinates.drop_duplicates(inplace=True)
    log_message(f"After dropping duplicates: {gene_coordinates.shape}")
    
    # Remove rows with NaN in start or end positions
    gene_coordinates.dropna(subset=[start_col, end_col], inplace=True)
    log_message(f"After removing NaN coordinates: {gene_coordinates.shape}")
    
    # Rename columns and create 'start' and 'end' columns
    gene_coordinates.rename(columns={chr_col: 'chr', 'ensgid': 'Gene'}, inplace=True)
    gene_coordinates['start'] = gene_coordinates[start_col].astype(int)
    gene_coordinates['end'] = gene_coordinates[end_col].astype(int)
    
    # Remove 'chr' prefix from the 'chr' column if present
    gene_coordinates['chr'] = gene_coordinates['chr'].str.replace('chr', '', regex=False)
    
    # Sort by 'chr' in natural chromosome order, then by 'start' and 'end'
    gene_coordinates['chr'] = pd.Categorical(
        gene_coordinates['chr'], 
        categories=[str(i) for i in range(1, 23)] + ['X', 'Y', 'M'],
        ordered=True
    )
    
    # Filter out chromosomes M, X, Y (keep only autosomes 1-22)
    gene_coordinates = gene_coordinates[~gene_coordinates['chr'].isin(['M', 'X', 'Y'])]
    gene_coordinates.sort_values(by=['chr', 'start', 'end'], inplace=True)
    log_message(f"After filtering autosomes 1-22: {gene_coordinates.shape}")
    
    # === Output 1: LDSC/mBAT version (chr, start, end, Gene) - NO HEADER ===
    log_message("Creating LDSC/mBAT gene coordinate file...")
    df_ldsc = gene_coordinates[['chr', 'start', 'end', 'Gene']].copy()
    df_ldsc.to_csv(args.output_ldsc, sep='\t', index=False, header=False)
    log_message(f"LDSC gene coordinate file saved to {args.output_ldsc} ({len(df_ldsc)} genes)")
    
    # === Output 2: MAGMA version (Gene, chr, start, end, strand, gene_name) - NO HEADER ===
    log_message("Creating MAGMA gene coordinate file...")
    df_magma = gene_coordinates[['Gene', 'chr', 'start', 'end', strand_col, 'gene_name']].copy()
    df_magma.to_csv(args.output_magma, sep='\t', index=False, header=False)
    log_message(f"MAGMA gene coordinate file saved to {args.output_magma} ({len(df_magma)} genes)")
    
    # === Output 3: CEPO/scDRS version (Gene, chr, start, end, gene_type, gene_name) - WITH HEADER ===
    log_message("Creating CEPO/scDRS gene coordinate file...")
    df_cepo = gene_coordinates[['Gene', 'chr', 'start', 'end', 'gene_type', 'gene_name']].copy()
    df_cepo.to_csv(args.output_cepo, sep='\t', index=False, header=True)
    log_message(f"CEPO gene coordinate file saved to {args.output_cepo} ({len(df_cepo)} genes)")
    
    log_message(f"Python script finished successfully for {args.genome_build}.")


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