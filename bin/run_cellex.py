#!/usr/bin/env python
"""
Compute CELLEX expression specificity metrics and extract the GES matrix
used by conLDSC-GES.
"""
import argparse
import time

import cellex
import numpy as np
import scanpy as sc
from scipy import sparse

GES_KEY = 'ges.esw_s'


def log_message(message):
    print(f"--- {time.strftime('%Y-%m-%d %H:%M:%S')} --- {message}", flush=True)


def assert_normalized(X):
    sample = X.data[:200000] if sparse.issparse(X) else np.asarray(X).ravel()[:200000]
    sample = sample[np.isfinite(sample)]
    if sample.size == 0:
        raise ValueError("Expression matrix appears to be empty.")
    if np.allclose(sample, np.round(sample)) and float(sample.max()) > 50.0:
        raise ValueError(
            "Expression matrix looks like raw integer counts, but CELLEX is being called "
            "with normalize=False. Run NORMALIZE_H5AD first so X holds log2(TPM+1)."
        )
    log_message(f"Expression scale check passed (max sampled value {float(sample.max()):.3f})")


def main(args):
    log_message(f"Reading {args.input_h5ad}")
    adata = sc.read_h5ad(args.input_h5ad)
    log_message(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes")

    if args.cell_type_col not in adata.obs.columns:
        raise ValueError(
            f"Cell type column {args.cell_type_col!r} not in obs. "
            f"Available: {sorted(adata.obs.columns)}"
        )

    assert_normalized(adata.X)

    log_message("Building CELLEX ESObject")
    eso = cellex.ESObject(
        data=adata.to_df().T,
        annotation=adata.obs[args.cell_type_col],
        normalize=False,
        verbose=True,
    )

    log_message("Computing CELLEX metrics")
    eso.compute(verbose=True)

    log_message(f"Saving all CELLEX metrics with prefix {args.prefix}")
    eso.save_as_csv(keys=["all"], path=".", file_prefix=args.prefix, verbose=True)

    if GES_KEY not in eso.results:
        raise KeyError(
            f"CELLEX did not produce {GES_KEY!r}. Available keys: {sorted(eso.results)}"
        )

    ges = eso.results[GES_KEY].copy()
    ges.index.name = 'gene'
    ges = ges.fillna(0.0)
    ges = ges.reset_index()

    ges.to_csv(args.output_ges, index=False)
    log_message(f"Wrote {args.output_ges}: {len(ges)} genes x {ges.shape[1] - 1} cell types")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input_h5ad', required=True)
    parser.add_argument('--cell_type_col', required=True)
    parser.add_argument('--prefix', required=True)
    parser.add_argument('--output_ges', required=True)
    main(parser.parse_args())
