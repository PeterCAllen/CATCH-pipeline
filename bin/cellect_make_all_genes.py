#!/usr/bin/env python
"""Port of CELLECT scripts/make_all_genes_snake.py to an argv interface."""
import argparse

import pandas as pd


def main(args):
    specificity_df = pd.read_csv(args.spec_matrix, index_col='gene')
    all_genes = pd.DataFrame(index=specificity_df.index)
    all_genes['all_genes_in_dataset'] = 1
    all_genes.to_csv(args.out)
    print(f"Wrote {args.out}: {len(all_genes)} genes")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--spec_matrix', required=True)
    parser.add_argument('--out', required=True)
    main(parser.parse_args())
