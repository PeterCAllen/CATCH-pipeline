# Legacy modules

Nothing in this directory is part of the CATCH combination as defined in
Li et al., *Benchmarking methods integrating GWAS and single-cell transcriptomic data
for mapping trait-cell type associations*. CATCH is the Cauchy combination of exactly
four components:

1. conLDSC-Cepo — `subworkflows/conldsc.nf`
2. conLDSC-GES — `subworkflows/conldsc.nf`
3. seismic-mBAT-combo — `subworkflows/seismic.nf`
4. scDRS — `subworkflows/scdrs.nf`

These modules are retained because they are useful standalone and because the MAGMA
path is needed to validate the seismic branch (running `seismicGWAS` against a MAGMA
`.genes.out` before substituting mBAT-combo z-statistics). They are **not** wired into
`COMBINE_CAUCHY` and are off by default.

## Contents

| File | What it was | Why it left the combination |
|---|---|---|
| `run_magma.nf`, `create_magma_geneset.nf`, `../../subworkflows/legacy/magma.nf` | MAGMA-GSEA. Was Cauchy input #2. | MAGMA-GSEA is benchmarked in the paper but is not one of the four CATCH components. |
| `create_continuous_beds.nf`, `create_ldsc_annot.nf`, `compute_ldsc_scores.nf`, `run_sldsc.nf`, `../../subworkflows/legacy/ldsc.nf` | The previous continuous-LDSC implementation: per-cell-type `ldsc.py --h2 --overlap-annot` on full annots, `GenomicRanges` overlaps, no all-genes control annotation. | Superseded by the CELLECT-LDSC port in `subworkflows/conldsc.nf`, which is the implementation the manuscript used. |
| `quantile_analysis.nf` | Gazal quantile machinery producing the top-quintile heritability enrichment p-value. | The paper states τ (coefficient) P-values "showed better trait-cell type prioritization than heritability enrichment P-values", so the enrichment statistic is not used. |
| `make_ldsc_annots.nf` | Binary LDSC annotations (`make_annot.py --thin-annot`). | Never wired in, references a missing `bin/make_annot.py`, and hardcodes hg38. Continuous annotations outperform binary per the paper. |
| `make_ldsc_continuous_annot.nf` | Fused annot + LD score process. | Never wired in. |
| `preprocess_h5ad.nf` | Donor-balanced cell subsampling. | Never wired in. Replaced by `modules/normalize_h5ad.nf`, which implements the manuscript's log2(TPM+1) contract. |

## Enabling MAGMA for seismic validation

```
nextflow run main.nf --run_magma true --run_conldsc false --run_seismic false --run_scdrs false ...
```

This produces `*.magma_genes.genes.out`, whose `GENE` and `ZSTAT` columns are the native
input to `seismicGWAS::get_ct_trait_associations`.
