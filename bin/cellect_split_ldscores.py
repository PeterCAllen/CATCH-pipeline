#!/usr/bin/env python
"""Port of CELLECT scripts/split_ldscores_snake.py to an argv interface."""
import argparse
import os
import re

import pandas as pd


def split_M_files(file_ldscore, chromosome, prefix_genomic_annot, list_annotations, out_dir):
    print(f"CHR={chromosome} | Writing .M and .M_5_50 files")
    file_ldscore_base = re.sub(r"\.l2\.ldscore\.gz$", "", file_ldscore)
    file_M = f"{file_ldscore_base}.l2.M"
    file_M_5_50 = f"{file_ldscore_base}.l2.M_5_50"
    with open(file_M, "r") as fh_M, open(file_M_5_50, "r") as fh_M_5_50:
        list_M = fh_M.readline().rstrip().split()
        list_M_5_50 = fh_M_5_50.readline().rstrip().split()
    assert len(list_M) == len(list_M_5_50) == len(list_annotations)
    for i in range(len(list_annotations)):
        annotation = list_annotations[i]
        file_out_M = f"{out_dir}/{prefix_genomic_annot}__{annotation}.{chromosome}.l2.M"
        file_out_M_5_50 = f"{out_dir}/{prefix_genomic_annot}__{annotation}.{chromosome}.l2.M_5_50"
        with open(file_out_M, "w") as fh_out_M, open(file_out_M_5_50, "w") as fh_out_M_5_50:
            fh_out_M.write(list_M[i])
            fh_out_M_5_50.write(list_M_5_50[i])
    print(f"CHR={chromosome} | DONE writing .M and .M_5_50 files")


def split_annot_file(file_ldscore, chromosome, prefix_genomic_annot, list_annotations, out_dir):
    print(f"CHR={chromosome} | START splitting annot file")
    file_ldscore_base = re.sub(r"\.l2\.ldscore\.gz$", "", file_ldscore)
    file_annot = f"{file_ldscore_base}.annot.gz"
    df_annot = pd.read_csv(file_annot, sep="\t")
    for i in range(len(list_annotations)):
        annotation = list_annotations[i]
        file_out_annot = f"{out_dir}/{prefix_genomic_annot}__{annotation}.{chromosome}.annot.gz"
        df_annot[[annotation]].to_csv(file_out_annot, sep="\t", index=False, compression="gzip")
    print(f"CHR={chromosome} | DONE splitting annot file")


def split_ldscore_file_per_annotation(file_ldscore, out_dir):
    print(f"Processing file_ldscore {file_ldscore}")
    m = re.search(r"(.*)\.COMBINED_ANNOT\.(\d{1,2})\.l2.ldscore.gz$", os.path.basename(file_ldscore))
    if m is None:
        raise ValueError(
            f"{file_ldscore} does not match the expected "
            "<prefix>.COMBINED_ANNOT.<chr>.l2.ldscore.gz naming."
        )
    prefix_genomic_annot = m.groups()[0]
    chromosome = m.groups()[1]

    df = pd.read_csv(file_ldscore, sep="\t")
    annotations_header = df.columns[3:].tolist()
    annotations_clean = [re.sub(r"L2$", "", x) for x in annotations_header]

    split_M_files(file_ldscore, chromosome, prefix_genomic_annot, annotations_clean, out_dir)
    split_annot_file(file_ldscore, chromosome, prefix_genomic_annot, annotations_clean, out_dir)

    for annotation in annotations_header:
        annotation_clean = re.sub(r"L2$", "", annotation)
        file_out_ldscore = (
            f"{out_dir}/{prefix_genomic_annot}__{annotation_clean}.{chromosome}.l2.ldscore.gz"
        )
        df[["CHR", "SNP", "BP", annotation]].to_csv(
            file_out_ldscore, sep="\t", index=False, compression="gzip"
        )

    print(f"CHR={chromosome} | Split {len(annotations_header)} annotations")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ldscore', required=True)
    parser.add_argument('--out_dir', default='.')
    a = parser.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    split_ldscore_file_per_annotation(a.ldscore, a.out_dir)
