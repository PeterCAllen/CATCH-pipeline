#!/usr/bin/env python
"""
Standardize a single-cell h5ad to the CATCH expression contract:
protein-coding autosomal genes, cell types with >= min_cells cells,
CELLECT-safe cell type labels, and log2(TPM + 1) values in X.
"""
import argparse
import re
import time

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

TARGET_SUM = 1e6
CELLECT_SAFE = re.compile(r'^[A-Za-z0-9_-]+$')


def log_message(message):
    print(f"--- {time.strftime('%Y-%m-%d %H:%M:%S')} --- {message}", flush=True)


def sanitize_label(label):
    s = re.sub(r'[^A-Za-z0-9_-]', '_', str(label))
    s = re.sub(r'_+', '_', s)
    s = s.strip('_-')
    return s


def sanitize_cell_types(adata, cell_type_col, mapping_out):
    original = adata.obs[cell_type_col].astype(str)
    uniq = sorted(original.unique())

    mapping = {}
    collisions = {}
    for label in uniq:
        clean = sanitize_label(label)
        if not clean:
            raise ValueError(f"Cell type label {label!r} sanitizes to an empty string.")
        if not CELLECT_SAFE.match(clean) or '__' in clean:
            raise ValueError(f"Cell type label {label!r} -> {clean!r} is not CELLECT-safe.")
        collisions.setdefault(clean, []).append(label)
        mapping[label] = clean

    clashed = {k: v for k, v in collisions.items() if len(v) > 1}
    if clashed:
        for clean, originals in clashed.items():
            log_message(f"ERROR: labels {originals} all sanitize to {clean!r}")
        raise ValueError(
            "Cell type label sanitization is not injective. Rename the conflicting labels "
            "in the input h5ad before running."
        )

    pd.DataFrame(
        {'original': list(mapping.keys()), 'sanitized': list(mapping.values())}
    ).to_csv(mapping_out, sep='\t', index=False)

    n_changed = sum(1 for k, v in mapping.items() if k != v)
    log_message(f"Sanitized {n_changed}/{len(mapping)} cell type labels; mapping written to {mapping_out}")

    adata.obs[cell_type_col] = pd.Categorical(original.map(mapping))
    return adata


def detect_scale(X):
    sample = X.data[:200000] if sparse.issparse(X) else np.asarray(X).ravel()[:200000]
    sample = sample[np.isfinite(sample)]
    if sample.size == 0:
        return 'raw_counts'
    is_integral = np.allclose(sample, np.round(sample))
    if is_integral:
        return 'raw_counts'
    if float(sample.max()) < 50.0:
        return 'lognorm'
    return 'raw_counts'


def apply_elementwise(X, fn):
    if sparse.issparse(X):
        X = X.tocsr(copy=True)
        X.data = fn(X.data)
        return X
    return fn(np.asarray(X, dtype=np.float64))


def main(args):
    log_message(f"Reading {args.input_h5ad}")
    adata = sc.read_h5ad(args.input_h5ad)
    log_message(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes")

    if args.use_raw:
        if adata.raw is None:
            raise ValueError("--use_raw requested but adata.raw is None")
        log_message(f"Promoting adata.raw ({adata.raw.shape[1]} genes) to X")
        adata = adata.raw.to_adata()

    if args.cell_type_col not in adata.obs.columns:
        raise ValueError(
            f"Cell type column {args.cell_type_col!r} not in obs. "
            f"Available: {sorted(adata.obs.columns)}"
        )

    adata.var_names_make_unique()

    if args.gene_coords:
        coords = pd.read_csv(args.gene_coords, sep='\t')
        keep_genes = set(coords['Gene'].astype(str))
        before = adata.n_vars
        mask = adata.var_names.astype(str).isin(keep_genes)
        n_match = int(mask.sum())
        if n_match == 0:
            raise ValueError(
                "No h5ad var_names matched the gene coordinate file. Expected Ensembl gene IDs "
                f"(e.g. {list(keep_genes)[:3]}) but saw {list(adata.var_names[:3])}."
            )
        adata = adata[:, mask].copy()
        log_message(f"Gene filter: {before} -> {adata.n_vars} genes matching the coordinate file")

    counts = adata.obs[args.cell_type_col].value_counts()
    keep_types = counts[counts >= args.min_cells].index
    dropped = counts[counts < args.min_cells]
    if len(dropped) > 0:
        log_message(f"Dropping {len(dropped)} cell types with < {args.min_cells} cells: {list(dropped.index)}")
    if len(keep_types) == 0:
        raise ValueError(f"No cell types have >= {args.min_cells} cells.")
    adata = adata[adata.obs[args.cell_type_col].isin(keep_types)].copy()
    adata.obs[args.cell_type_col] = adata.obs[args.cell_type_col].astype(str)
    log_message(f"Cell type filter: {len(keep_types)} cell types, {adata.n_obs} cells retained")

    adata = sanitize_cell_types(adata, args.cell_type_col, args.output_mapping)

    scale = args.input_scale
    if scale == 'auto':
        scale = detect_scale(adata.X)
        log_message(f"Detected input scale: {scale}")
    else:
        log_message(f"Input scale (declared): {scale}")

    if scale == 'lognorm':
        if args.log_base == '2':
            log_message("Reversing prior log2(x+1) normalization")
            adata.X = apply_elementwise(adata.X, lambda d: np.exp2(d) - 1.0)
        else:
            log_message("Reversing prior log1p (natural) normalization")
            adata.X = apply_elementwise(adata.X, np.expm1)

    log_message(f"Normalizing each cell to {TARGET_SUM:.0e} total counts")
    sc.pp.normalize_total(adata, target_sum=TARGET_SUM)

    log_message("Applying log2(TPM + 1)")
    adata.X = apply_elementwise(adata.X, lambda d: np.log2(d + 1.0))

    adata.uns['catch_normalization'] = {
        'scheme': 'log2(TPM+1)',
        'target_sum': TARGET_SUM,
        'input_scale': scale,
        'min_cells': args.min_cells,
        'cell_type_col': args.cell_type_col,
    }

    log_message(f"Writing {args.output_h5ad}")
    adata.write_h5ad(args.output_h5ad, compression='gzip')

    ct_counts = adata.obs[args.cell_type_col].value_counts().sort_index()
    ct_counts.rename_axis('cell_type').reset_index(name='n_cells').to_csv(
        args.output_cell_counts, sep='\t', index=False
    )

    log_message(f"Final: {adata.n_obs} cells x {adata.n_vars} genes, {len(ct_counts)} cell types")
    log_message("Done.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input_h5ad', required=True)
    parser.add_argument('--output_h5ad', required=True)
    parser.add_argument('--output_mapping', required=True)
    parser.add_argument('--output_cell_counts', required=True)
    parser.add_argument('--cell_type_col', required=True)
    parser.add_argument('--gene_coords', default=None)
    parser.add_argument('--min_cells', type=int, default=20)
    parser.add_argument('--input_scale', choices=['auto', 'raw_counts', 'lognorm'], default='auto')
    parser.add_argument('--log_base', choices=['e', '2'], default='e')
    parser.add_argument('--use_raw', action='store_true')
    main(parser.parse_args())
