# CATCH Pipeline

Nextflow implementation of CATCH (Li et al.), which combines two methods via Cauchy
combination to map trait–cell type associations from GWAS and scRNA-seq data:

1. **conLDSC-Cepo** — CELLECT-LDSC on a Cepo specificity matrix
2. **scDRS** — scored on an mBAT-combo gene set

Two further methods remain available as optional, off-by-default validation branches
(not part of the CATCH combination itself):

- **conLDSC-GES** — CELLECT-LDSC on the CELLEX GES matrix (`--run_cellex true`)
- **seismic-mBAT-combo** — `seismicGWAS` using mBAT-combo z-statistics (`--run_seismic true`)

MAGMA-GSEA and binary LDSC are in `modules/legacy/` and are off by default.

## Requirements

- Nextflow >= 21.04.0
- Singularity
- A node with ~256 GB RAM (Cepo and CELLEX)
- hg19/GRCh37 reference data (1000G EUR Phase 3, 53-annotation baseline, CELLECT's
  `print_snps.txt`, `GRCh37-chr-sizes.txt`, `w_hm3.snplist`) under `data/reference/` --
  freely downloadable from the [S-LDSC reference files Zenodo record](https://doi.org/10.5281/zenodo.7768714)
  (`1000G_Phase3_plinkfiles.tgz`, `1000G_Phase3_weights_hm3_no_MHC.tgz`,
  `1000G_Phase3_baseline_v1.2_ldscores.tgz` -- note: v1.2, not the original v1.1_thin_annot,
  which Broad Institute moved to a requester-pays bucket) and from
  [perslab/CELLECT](https://github.com/perslab/CELLECT) (`data/ldsc/{print_snps.txt,GRCh37-chr-sizes.txt,w_hm3.snplist}`,
  served via Git LFS at `media.githubusercontent.com`). `data/` is gitignored -- keep a
  copy on non-purged storage (e.g. `/g/data` on NCI systems), not `/scratch`.

Build the containers first (requires root or `--fakeroot`/subuid access -- on clusters
without that, build on your own machine with Docker/root and copy the resulting `.sif`
files onto a non-purged path such as `/g/data`, not `/scratch`):

```bash
cd environments
for d in ldsc-py3 cellect-py3 cellex seismic py-r-cepo-scdrs gcta_v1.94.1 scdrs_v1.0.2; do
    singularity build --fakeroot $d.sif $d.def
done
```

`py-r-cepo-scdrs.def`, `gcta_v1.94.1.def`, and `scdrs_v1.0.2.def` were authored by
inspecting how the pipeline calls into each image (see the `%labels` in each file);
`gcta_v1.94.1.def` and `scdrs_v1.0.2.def` pin exact upstream versions, but
`py-r-cepo-scdrs.def`'s R/Python package versions are best-effort and were not
verified against a prior working image -- check results against any known-good
reference before trusting them for publication numbers.

## Usage

```bash
nextflow run main.nf \
    --h5ad_input data/sc.h5ad \
    --gwas_sumstats data/gwas.txt.gz \
    --gwas_name AD_Jansen2019 \
    --genome_build hg19 \
    --gwas_sample_size 455258 \
    --gene_matrix data/geneMatrix.tsv.gz \
    --cell_type_col cell.labels \
    --outdir results/AD_Jansen2019 \
    -profile pbs
```

`nextflow run main.nf --help` lists all parameters.

If you already have pre-formatted GWAS files (as consortium GWAS are often
distributed), skip the pipeline's own reformatting and pass them directly instead of
`--gwas_sumstats`:
- `--gwas_cojo <trait>.ma` — GCTA-COJO format (`SNP A1 A2 freq b se p N`), bypasses
  reformatting for the seismic/scDRS (mBAT-combo) branch.
- `--gwas_sumstats_munged <trait>.sumstats.gz` — pre-munged LDSC format, bypasses
  `munge_sumstats.py` for the conLDSC branch.

At least one of `--gwas_sumstats` / `--gwas_cojo` / `--gwas_sumstats_munged` is required,
depending on which components you run.

Components can be run individually, for example `--run_cellex true --run_scdrs false`.
The Cauchy combination (CATCH) runs only when `--run_conldsc --run_cepo --run_scdrs`
are all enabled; `--run_cellex`/`--run_seismic` are optional validation branches and
are not required for (or included in) the combination.

The paper applies FDR across all cell types *and traits*. One run covers one trait, so
after running every trait:

```bash
Rscript bin/catch_fdr.R results/catch_fdr.tsv results/*/combined/*_catch_combined.tsv
```

## Output

```
<outdir>/
├── preprocessed/   normalized h5ad, cell type mapping
├── metrics/        Cepo (and, if enabled, CELLEX) specificity matrices
├── conldsc/        prioritization.csv per specificity matrix
├── mbat/           mBAT-combo results
├── seismic/        <gwas>_seismic.tsv (only if --run_seismic true)
├── scdrs/          cell scores and group results
├── combined/       <gwas>_catch_combined.tsv
└── pipeline_info/  execution_trace/report/timeline (compute-hours, peak memory, etc.)
```

`combined/<gwas>_catch_combined.tsv` is the main result: the two component p-values and
`CATCH_P` for each cell type.

## Layout

| Path | Purpose |
|---|---|
| `subworkflows/metrics.nf` | Cepo and CELLEX/GES |
| `subworkflows/conldsc.nf` | CELLECT-LDSC port |
| `subworkflows/mbat.nf` | GCTA mBAT-combo, shared by seismic and scDRS |
| `subworkflows/seismic.nf` | seismic-mBAT-combo |
| `subworkflows/scdrs.nf` | scDRS |
| `modules/legacy/` | MAGMA-GSEA and the superseded LDSC path |
| `tests/` | tests for the CELLECT ports |

## Citation

Cite the CATCH paper along with the component methods: LDSC (Bulik-Sullivan 2015),
CELLECT and CELLEX (Timshel 2020), Cepo (Kim 2021), seismic (Lai 2025),
mBAT-combo (Li 2023), scDRS (Zhang 2022).

## Authors

Ang Li, Jian Zeng, Peter C Allen

Issues: https://github.com/PeterCAllen/CATCH-pipeline/issues
