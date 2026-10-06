#!/usr/bin/env python
"""
Validate and normalize a gene x cell-type specificity matrix into the form
CELLECT expects: comma-separated, first column named 'gene', no NAs,
no all-zero annotations, and CELLECT-safe annotation names.
"""
import argparse
import re
import time

import pandas as pd

CELLECT_SAFE = re.compile(r'^[A-Za-z0-9_-]+$')


def log_message(message):
    print(f"--- {time.strftime('%Y-%m-%d %H:%M:%S')} --- {message}", flush=True)


def main(args):
    log_message(f"Reading {args.input}")
    sep = '\t' if args.input_sep == 'tab' else ','
    df = pd.read_csv(args.input, sep=sep)

    first = df.columns[0]
    if first != 'gene':
        log_message(f"Renaming first column {first!r} -> 'gene'")
        df = df.rename(columns={first: 'gene'})

    df['gene'] = df['gene'].astype(str).str.strip()
    df = df[df['gene'].ne('') & df['gene'].ne('nan')]

    n_dup = int(df['gene'].duplicated().sum())
    if n_dup:
        log_message(f"WARNING: {n_dup} duplicate gene IDs; keeping the first occurrence of each")
        df = df[~df['gene'].duplicated(keep='first')]

    annots = [c for c in df.columns if c != 'gene']
    if not annots:
        raise ValueError("No annotation columns found after the gene column.")

    bad = [c for c in annots if not CELLECT_SAFE.match(str(c)) or '__' in str(c)]
    if bad:
        raise ValueError(
            f"Annotation names are not CELLECT-safe: {bad}. "
            "Cell type labels should have been sanitized by normalize_h5ad.py; "
            "check that the metric ran on the normalized h5ad."
        )

    df[annots] = df[annots].apply(pd.to_numeric, errors='coerce')

    n_na = int(df[annots].isna().sum().sum())
    if n_na:
        log_message(f"Filling {n_na} NA values with 0")
        df[annots] = df[annots].fillna(0.0)

    zero_cols = [c for c in annots if float(df[c].abs().sum()) == 0.0]
    if zero_cols:
        log_message(f"Dropping {len(zero_cols)} all-zero annotations: {zero_cols}")
        df = df.drop(columns=zero_cols)
        annots = [c for c in annots if c not in zero_cols]

    if not annots:
        raise ValueError("Every annotation column was all-zero; nothing left to analyze.")

    if args.top_pct is not None:
        if not 0 < args.top_pct <= 1:
            raise ValueError(f"--top_pct must be in (0, 1], got {args.top_pct}")
        n_keep = max(1, int(round(len(df) * args.top_pct)))
        log_message(f"Restricting each annotation to its top {args.top_pct:.0%} of genes (n={n_keep})")
        for c in annots:
            cutoff = df[c].nlargest(n_keep).min()
            df.loc[df[c] < cutoff, c] = 0.0

    zero_genes = int((df[annots].abs().sum(axis=1) == 0).sum())
    log_message(f"{zero_genes}/{len(df)} genes are zero across all annotations (retained; CELLECT tolerates these)")

    df = df[['gene'] + annots]
    df.to_csv(args.output, index=False)

    log_message(f"Wrote {args.output}: {len(df)} genes x {len(annots)} annotations")
    with open(args.output_annotations, 'w') as fh:
        fh.write('\n'.join(annots) + '\n')
    log_message(f"Wrote annotation list to {args.output_annotations}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--output_annotations', required=True)
    parser.add_argument('--input_sep', choices=['comma', 'tab'], default='comma')
    parser.add_argument('--top_pct', type=float, default=None)
    main(parser.parse_args())
