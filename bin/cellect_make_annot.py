#!/usr/bin/env python
"""
Port of CELLECT scripts/make_annot_from_geneset_all_chr_snake.py to an argv interface.

For continuous annotations, a variant spanned by multiple genes within the window is
assigned the maximum annotation value.
"""
import argparse
import os

import numpy as np
import pandas as pd
import pybedtools


def get_max_ES_vals(overlap_df, cell_type_ES_dict):
    max_vals = []
    for segment in overlap_df['GENES']:
        gene_list = segment.split(';')
        ES_val_list = []
        for gene in gene_list:
            if gene in cell_type_ES_dict.keys():
                ES_val_list.append(cell_type_ES_dict[gene])
            else:
                ES_val_list.append(0)
        max_ES_val = np.max(ES_val_list)
        max_vals.append(max_ES_val)
    segment_df = overlap_df.copy()
    segment_df['SCORE'] = max_vals
    segment_df['CHR'] = 'chr' + segment_df['CHR'].astype(str)
    return segment_df


def map_ES_values_to_genes(chr_overlap_df, specificity_df):
    ES_dict = {}
    for cell_type in specificity_df.columns:
        ES_dict[cell_type] = specificity_df[cell_type].to_dict()

    max_ES_vals_dict = {}
    for cell_type in specificity_df.columns:
        cell_type_ES_dict = ES_dict[cell_type]
        max_ES_vals_dict[cell_type] = pybedtools.bedtool.BedTool.from_dataframe(
            get_max_ES_vals(chr_overlap_df, cell_type_ES_dict)
        )
    return max_ES_vals_dict


def make_annot_file_per_chromosome(chromosome, dict_of_beds, out_dir, out_prefix,
                                   bimfile, all_genes, keep_annots, out_annot_keep):
    print(f'making annot files for chromosome {chromosome}')
    df_bim = pd.read_csv(bimfile, sep=r'\s+', usecols=[0, 1, 2, 3],
                         names=['CHR', 'SNP', 'CM', 'BP'])

    iter_bim = [['chr' + str(x1), str(x2), str(x2)] for (x1, x2) in np.array(df_bim[['CHR', 'BP']])]
    bimbed = pybedtools.BedTool(iter_bim)

    counter = 1
    list_df_annot = []
    annotation_genes_df = None

    for name_annotation in sorted(dict_of_beds):
        print(f"CHR={chromosome} | annotation={name_annotation}, #{counter}/#{len(dict_of_beds)}")
        bed_for_annot = dict_of_beds[name_annotation]
        annotbed = bimbed.intersect(bed_for_annot, wb=True)

        bp = [x.start for x in annotbed]
        annotation_value = [x.fields[7] for x in annotbed]

        if keep_annots and counter == 1:
            annotation_genes = [x.fields[6] for x in annotbed]
            annotation_genes_df = pd.DataFrame(
                {'BP': bp, 'GENES': annotation_genes}
            ).drop_duplicates(subset='BP', keep='first')

        os.remove(annotbed.fn)

        df_annot_overlap_bp = pd.DataFrame({'BP': bp, name_annotation: annotation_value})
        df_annot_overlap_bp = df_annot_overlap_bp.drop_duplicates(subset='BP', keep='first')
        df_annot = pd.merge(df_bim, df_annot_overlap_bp, how='left', on='BP')
        df_annot = df_annot[[name_annotation]]
        df_annot.fillna(0, inplace=True)
        df_annot[name_annotation] = df_annot[name_annotation].astype(float)
        list_df_annot.append(df_annot)
        counter += 1

    print(f"CHR={chromosome} | Concatenating annotations...")
    df_annot_combined = pd.concat(list_df_annot, axis='columns')

    print(f"CHR={chromosome} | Writing annotations...")
    if all_genes:
        file_out_annot_combined = f"{out_dir}/all_genes_in_{out_prefix}.{chromosome}.annot.gz"
    else:
        file_out_annot_combined = f"{out_dir}/{out_prefix}.COMBINED_ANNOT.{chromosome}.annot.gz"

    if keep_annots:
        df_annot_keep_combined = pd.concat([df_bim] + list_df_annot, axis='columns')
        df_annot_keep_combined = pd.merge(df_annot_keep_combined, annotation_genes_df,
                                          how='left', on='BP')
        df_annot_keep_combined.to_csv(out_annot_keep, sep="\t", index=False, compression="gzip")

    df_annot_combined.to_csv(file_out_annot_combined, sep="\t", index=False, compression="gzip")
    print(f"Wrote {file_out_annot_combined}")


def main(args):
    overlap_df = pd.read_csv(args.overlap_segments, sep=r'\s+',
                             names=['CHR', 'START', 'END', 'GENES'])
    specificity_df = pd.read_csv(args.spec_matrix, index_col='gene')

    if specificity_df.isna().any().any():
        raise ValueError(
            "Specificity matrix contains NA values; np.max would propagate NaN into the "
            "annotation. Run sanitize_specificity_matrix.py first."
        )

    dict_of_beds = map_ES_values_to_genes(overlap_df, specificity_df)

    make_annot_file_per_chromosome(
        args.chromosome, dict_of_beds, args.out_dir, args.out_prefix,
        args.bimfile, args.all_genes, args.keep_annots, args.out_annot_keep
    )
    print("Make annot script is done!")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--spec_matrix', required=True)
    parser.add_argument('--overlap_segments', required=True)
    parser.add_argument('--bimfile', required=True)
    parser.add_argument('--chromosome', required=True)
    parser.add_argument('--out_dir', default='.')
    parser.add_argument('--out_prefix', required=True)
    parser.add_argument('--all_genes', action='store_true')
    parser.add_argument('--keep_annots', action='store_true')
    parser.add_argument('--out_annot_keep', default=None)
    main(parser.parse_args())
