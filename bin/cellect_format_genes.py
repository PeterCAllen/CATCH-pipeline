#!/usr/bin/env python
"""Port of CELLECT scripts/format_and_map_snake.py to an argv interface."""
import argparse

import numpy as np
import pandas as pd


def format_genes(gene_coord_df, chr_sizes_df, out_dir, windowsize_kb):
    chr_sizes_dict = chr_sizes_df.to_dict()[1]
    windowsize = windowsize_kb * 1000

    gene_coord_df = gene_coord_df[gene_coord_df['CHR'].isin(chr_sizes_dict.keys())]
    gene_coord_df = gene_coord_df.dropna().sort_values(by=['CHR', 'START'])

    gene_coord_df['START'] = np.maximum(0, gene_coord_df['START'] - windowsize)
    gene_coord_df['END'] = np.minimum(
        gene_coord_df['END'] + windowsize, gene_coord_df['CHR'].map(chr_sizes_dict)
    )

    for chrom_num in chr_sizes_dict.keys():
        chrom_df = gene_coord_df[gene_coord_df['CHR'] == chrom_num]
        chrom_df = chrom_df[['CHR', 'START', 'END', 'GENE']]
        out_file = f"{out_dir}/genes_plus_{windowsize_kb}kb.{chrom_num}.bed"
        chrom_df.to_csv(out_file, sep='\t', index=False, header=False)
        print(f"chr{chrom_num}: {len(chrom_df)} genes -> {out_file}")


def main(args):
    gene_coord_df = pd.read_csv(
        args.gene_coords, sep=r'\s+', index_col=None, header=None,
        names=['GENE', 'CHR', 'START', 'END', 'STRAND', 'GENE_NAME']
    )
    gene_coord_df = gene_coord_df.drop(['STRAND', 'GENE_NAME'], axis=1)

    chr_sizes_df = pd.read_csv(args.chr_sizes, sep=r'\s+', index_col=0, header=None)

    format_genes(gene_coord_df, chr_sizes_df, args.out_dir, args.windowsize_kb)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gene_coords', required=True)
    parser.add_argument('--chr_sizes', required=True)
    parser.add_argument('--out_dir', required=True)
    parser.add_argument('--windowsize_kb', type=int, required=True)
    main(parser.parse_args())
