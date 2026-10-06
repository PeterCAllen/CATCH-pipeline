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
    
    # --- Take 10k max for each cell-type ---
    np.random.seed(602)

    cells_per_type = 10000
    celltype_col = args.cell_type_col
    donor_col = 'individual'  

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
    
    args = parser.parse_args()
    main(args)
