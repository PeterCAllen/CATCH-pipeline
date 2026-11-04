#!/usr/bin/env python
import scanpy as sc
import numpy as np
import argparse

def main(args):
    """
    Subsets, filters, and normalizes a raw h5ad file using backed mode for memory efficiency.
    """
    print("Reading h5ad file in backed mode (memory-efficient)...")
    # Use backed mode to avoid loading entire dataset into memory
    adata = sc.read_h5ad(args.input_h5ad, backed='r')
    print(f"Dimensions before filtering: {adata.n_obs} cells x {adata.n_vars} genes")

    # Get metadata for sampling without loading the full matrix
    obs_df = adata.obs.copy()

    # --- Take at most cells_per_type cells for each cell-type ---
    np.random.seed(602)

    cells_per_type = args.cells_per_type
    celltype_col = args.cell_type_col
    donor_col = args.donor_col

    if args.cell_types:
        missing = set(args.cell_types) - set(obs_df[celltype_col].unique())
        if missing:
            raise ValueError(f"Cell types not found in '{celltype_col}': {sorted(missing)}")
        obs_df = obs_df[obs_df[celltype_col].isin(args.cell_types)]

    sampled_indices = []

    for cell_type in obs_df[celltype_col].unique():
        celltype_df = obs_df[obs_df[celltype_col] == cell_type]
        donors = celltype_df[donor_col].unique()
        n_per_donor = cells_per_type // len(donors)
        donor_sampled = []
        for donor in donors:
            donor_idx = celltype_df[celltype_df[donor_col] == donor].index
            n = min(n_per_donor, len(donor_idx))
            sampled = np.random.choice(donor_idx, n, replace=False)
            donor_sampled.extend(sampled)
        # If not enough cells (due to rounding), sample remaining randomly from all
        if len(donor_sampled) < cells_per_type:
            remaining = list(set(celltype_df.index) - set(donor_sampled))
            n_remaining = min(cells_per_type - len(donor_sampled), len(remaining))
            if n_remaining > 0:
                donor_sampled.extend(np.random.choice(remaining, n_remaining, replace=False))
        sampled_indices.extend(donor_sampled)

    print(f"Selected {len(sampled_indices)} cells for subsetting")
    
    # Subset using backed mode - this loads only the selected data
    print("Loading selected cells...")
    ad_sub = adata[sampled_indices].to_memory()
    
    # Convert raw to main if available
    if ad_sub.raw is not None:
        print("Converting raw data to main matrix...")
        ad_main = ad_sub.raw.to_adata()
        ad_main.layers['counts'] = ad_main.X.copy() # for R compatibility
        ad_sub = ad_main
    else:
        ad_sub.layers['counts'] = ad_sub.X.copy()
    
    # The pipeline matches genes by Ensembl ID, so var_names must be Ensembl IDs.
    if args.gene_id_col:
        if args.gene_id_col not in ad_sub.var.columns:
            raise ValueError(f"'{args.gene_id_col}' not in var. Available: {list(ad_sub.var.columns)}")
        ad_sub.var['gene_symbol'] = ad_sub.var_names.astype(str)
        ad_sub.var_names = ad_sub.var[args.gene_id_col].astype(str).values
        ad_sub.var_names_make_unique()
        print(f"Set var_names from var['{args.gene_id_col}'] (symbols kept in var['gene_symbol'])")

    print(f"Dimensions after filtering: {ad_sub.n_obs} cells x {ad_sub.n_vars} genes")

    # --- Save output ---
    print("Saving filtered h5ad file...")
    ad_sub.write(args.output_h5ad)
    print(f"Filtered h5ad file saved to {args.output_h5ad}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Filter and normalize a raw h5ad file.")
    parser.add_argument('--input_h5ad', type=str, required=True, help='Path to the input raw h5ad file.')
    parser.add_argument('--output_h5ad', type=str, required=True, help='Path for the output processed h5ad file.')
    parser.add_argument('--cell_type_col', type=str, required=True, help='Metadata column containing cell type labels.')
    parser.add_argument('--donor_col', type=str, default='individual', help='Metadata column identifying the donor; cells are sampled evenly across donors.')
    parser.add_argument('--cells_per_type', type=int, default=10000, help='Maximum number of cells to keep per cell type.')
    parser.add_argument('--cell_types', type=str, nargs='+', default=None, help='Only keep these cell types (default: all).')
    parser.add_argument('--gene_id_col', type=str, default=None, help='var column holding Ensembl gene IDs to use as var_names (the pipeline requires Ensembl IDs).')

    args = parser.parse_args()
    main(args)
