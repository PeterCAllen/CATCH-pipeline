#!/usr/bin/env python
"""
Port of the prioritization block of CELLECT scripts/parse_results.py.

Converts an ldsc.py --h2-cts .cell_type_results.txt into CELLECT's
prioritization.csv schema: gwas, specificity_id, annotation, beta, beta_se, pvalue.
"""
import argparse

import pandas as pd


def main(args):
    prioritization = pd.read_csv(args.cell_type_results, sep='\t', header=0)

    expected = {'Name', 'Coefficient', 'Coefficient_std_error', 'Coefficient_P_value'}
    missing = expected - set(prioritization.columns)
    if missing:
        raise ValueError(
            f"{args.cell_type_results} is missing expected columns {sorted(missing)}; "
            f"found {list(prioritization.columns)}"
        )

    prioritization['gwas'] = args.gwas

    specificity_id_annotation = prioritization['Name'].str.split('__', expand=True)
    if specificity_id_annotation.shape[1] < 2:
        raise ValueError(
            "Could not split 'Name' on '__' into specificity_id and annotation. "
            "Cell type labels must not contain '__'."
        )
    prioritization['specificity_id'] = specificity_id_annotation[0]
    prioritization['annotation'] = specificity_id_annotation[1]

    prioritization = prioritization[
        ['gwas', 'specificity_id', 'annotation',
         'Coefficient', 'Coefficient_std_error', 'Coefficient_P_value']
    ]
    prioritization = prioritization.rename(columns={
        'Coefficient': 'beta',
        'Coefficient_std_error': 'beta_se',
        'Coefficient_P_value': 'pvalue',
    })
    prioritization = prioritization.sort_values(['gwas', 'specificity_id', 'pvalue'])

    prioritization.to_csv(args.out, index=False)
    print(f"Wrote {args.out}: {len(prioritization)} cell types")

    n_missing = int(prioritization['pvalue'].isna().sum())
    if n_missing:
        print(
            f"WARNING: {n_missing} cell types have no p-value. The LDSC fork logs these as "
            "'*CTS ERROR*' and omits them, usually from a zero-variance annotation."
        )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cell_type_results', required=True)
    parser.add_argument('--gwas', required=True)
    parser.add_argument('--out', required=True)
    main(parser.parse_args())
