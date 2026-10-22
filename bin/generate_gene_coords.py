#!/usr/bin/env python
"""
Generate genome-build-specific gene coordinate files from geneMatrix.tsv.gz
Based on manual analysis scripts: 0-Gene_coordinate_ready-hg19.py and 0-Gene_coordinate_ready-hg38.py
"""
import pandas as pd
import argparse
import shutil
import time

EXPECTED_PC_GENES = 19430


def log_message(message):
    """Prints a formatted log message with a timestamp."""
    timestamp = time.strftime('%Y-%m-%d %H:%M:%S')
    print(f"--- {timestamp} --- {message}", flush=True)

def protein_coding_mask(df):
    """Boolean mask of rows whose gene_type is protein_coding (case/dash-insensitive)."""
    if 'gene_type' not in df.columns:
        raise ValueError(
            "Protein-coding filter requested but the gene matrix has no 'gene_type' column. "
            f"Available columns: {sorted(df.columns)}"
        )
    gene_type_norm = df['gene_type'].astype(str).str.strip().str.lower().str.replace('-', '_', regex=False)
    return gene_type_norm == 'protein_coding'


def write_cellect_ref(cellect_ref, output, gene_matrix=None):
    """Pass the CELLECT gene coordinate reference through; if gene_matrix is given, keep only
    genes whose Ensembl ID is protein-coding in it."""
    if gene_matrix is None:
        shutil.copyfile(cellect_ref, output)
        log_message(f"CELLECT reference gene coordinates copied unchanged (all genes) to {output}")
        return

    ref = pd.read_csv(cellect_ref, sep=r'\s+', header=None, dtype=str)
    if ref.shape[1] != 6:
        raise ValueError(f"{cellect_ref}: expected 6 columns (GENE CHR START END STRAND GENE_NAME), got {ref.shape[1]}")
    in_pc = ref[0].isin(gene_matrix.loc[protein_coding_mask(gene_matrix), 'ensgid'])
    n_unmatched = (~ref[0].isin(gene_matrix['ensgid'])).sum()
    ref[in_pc].to_csv(output, sep='\t', index=False, header=False)
    log_message(
        f"CELLECT reference gene coordinates restricted to protein_coding: kept {in_pc.sum()}/{len(ref)} "
        f"(dropped {len(ref) - in_pc.sum() - n_unmatched} non-protein-coding, {n_unmatched} absent "
        f"from the gene matrix). Biotypes come from the gene matrix annotation, not the reference's own."
    )
    log_message(f"CELLECT reference gene coordinate file saved to {output}")


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

    if args.cellect_ref:
        write_cellect_ref(args.cellect_ref, args.output_cellect_ref,
                          gene_coordinates if args.coords_protein_coding else None)
    
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
    log_message(f"After restricting to autosomes 1-22: {df_coords.shape}")

    # mBAT-combo is a pure GWAS/LD-based gene test with no expression-quantification
    # rationale for excluding non-coding loci (unlike Cepo/scDRS, where protein-coding
    # restriction is about scRNA-seq measurement reliability). By default give mBAT the
    # full gene set, independent of --protein_coding/--no_protein_coding below;
    # --coords_protein_coding restricts it to protein-coding genes to save compute
    # (downstream scDRS/seismic only ever use genes present in the PC-filtered h5ad).
    df_mbat = df_coords.copy()
    if args.coords_protein_coding:
        df_mbat = df_mbat[protein_coding_mask(df_mbat)]
        log_message(f"mBAT gene coordinates (protein_coding only, --coords_protein_coding): {df_mbat.shape}")
    else:
        log_message(f"mBAT gene coordinates (unfiltered, all gene types): {df_mbat.shape}")
    df_mbat.sort_values(by=['chr', 'start', 'end'], inplace=True)

    if args.protein_coding:
        df_coords = df_coords[protein_coding_mask(df_coords)]
        log_message(f"After restricting to protein_coding: {df_coords.shape}")

        n_genes = len(df_coords)
        if not (EXPECTED_PC_GENES * 0.8) <= n_genes <= (EXPECTED_PC_GENES * 1.2):
            log_message(
                f"WARNING: {n_genes} protein-coding autosomal genes is >20% away from the "
                f"expected {EXPECTED_PC_GENES}. Check the gene matrix build and annotation version."
            )
        else:
            log_message(f"Gene universe size {n_genes} consistent with expected ~{EXPECTED_PC_GENES}.")
    else:
        log_message("Protein-coding filter DISABLED (--no_protein_coding).")

    df_coords.sort_values(by=['chr', 'start', 'end'], inplace=True)
    log_message(f"Gene coordinates (hg19): {df_coords.shape}")

    # === Output 1: LDSC version (chr, start, end, Gene) - NO HEADER - hg19 ===
    log_message("Creating LDSC gene coordinate file (hg19)...")
    df_ldsc_out = df_coords[['chr', 'start', 'end', 'Gene']].copy()
    df_ldsc_out.to_csv(args.output_ldsc, sep='\t', index=False, header=False)
    log_message(f"LDSC gene coordinate file saved to {args.output_ldsc} ({len(df_ldsc_out)} genes, hg19 coordinates)")

    # === Output mBAT: all genes unless --coords_protein_coding (chr, start, end, Gene) - NO HEADER ===
    log_message("Creating mBAT gene coordinate file (hg19)...")
    df_mbat_out = df_mbat[['chr', 'start', 'end', 'Gene']].copy()
    df_mbat_out.to_csv(args.output_mbat, sep='\t', index=False, header=False)
    log_message(f"mBAT gene coordinate file saved to {args.output_mbat} ({len(df_mbat_out)} genes, hg19 coordinates)")

    # === Output 2: MAGMA version (Gene, chr, start, end, strand, gene_name) - NO HEADER - hg19 ===
    log_message("Creating MAGMA gene coordinate file (hg19)...")
    df_magma_out = df_coords[['Gene', 'chr', 'start', 'end', coord_strand_col, 'gene_name']].copy()
    df_magma_out.to_csv(args.output_magma, sep='\t', index=False, header=False)
    log_message(f"MAGMA gene coordinate file saved to {args.output_magma} ({len(df_magma_out)} genes, hg19 coordinates)")
    
    # === Output 2b: CELLECT version (GENE, CHR, START, END, STRAND, GENE_NAME) - NO HEADER - hg19 ===
    if args.output_cellect:
        log_message("Creating CELLECT gene coordinate file (hg19)...")
        df_cellect_out = df_coords[['Gene', 'chr', 'start', 'end', coord_strand_col, 'gene_name']].copy()
        if args.coords_protein_coding and not args.protein_coding:
            df_cellect_out = df_cellect_out[protein_coding_mask(df_coords)]
        df_cellect_out.to_csv(args.output_cellect, sep='\t', index=False, header=False)
        log_message(f"CELLECT gene coordinate file saved to {args.output_cellect} ({len(df_cellect_out)} genes)")

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
    parser.add_argument('--output_mbat', type=str, required=True,
                        help='Path for the output mBAT gene location file (no header; all genes '
                             'unless --coords_protein_coding)')
    parser.add_argument('--output_cepo', type=str, required=True,
                        help='Path for the output CEPO gene coordinate file (with header)')
    parser.add_argument('--output_cellect', type=str, default=None,
                        help='Path for the output CELLECT gene coordinate file (no header)')
    parser.add_argument('--no_protein_coding', dest='protein_coding', action='store_false',
                        help='Disable the protein-coding gene filter')
    parser.set_defaults(protein_coding=True)
    parser.add_argument('--coords_protein_coding', action='store_true',
                        help='Restrict the mBAT, CELLECT and --cellect_ref gene coordinate files to '
                             'protein-coding genes (independent of --no_protein_coding)')
    parser.add_argument('--cellect_ref', type=str, default=None,
                        help='CELLECT-format gene coordinate reference (6 cols, no header) to pass '
                             'through to --output_cellect_ref, PC-filtered if --coords_protein_coding')
    parser.add_argument('--output_cellect_ref', type=str, default=None,
                        help='Path for the (optionally PC-filtered) copy of --cellect_ref')

    args = parser.parse_args()
    if bool(args.cellect_ref) != bool(args.output_cellect_ref):
        parser.error('--cellect_ref and --output_cellect_ref must be given together')
    main(args)