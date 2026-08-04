# CATCH Pipeline

Nextflow implementation of **CATCH** — *Cell-type Associations with Traits via CaucHy
combination of strategies* — from Li et al., *Benchmarking methods integrating GWAS and
single-cell transcriptomic data for mapping trait-cell type associations*.

CATCH is the Cauchy combination of exactly four component methods, selected in the paper
from 184,509 candidate combinations because they span all three methodological families:

| # | Component | Family | Implementation |
|---|---|---|---|
| 1 | **conLDSC-Cepo** | SC-to-GWAS | CELLECT-LDSC on a Cepo specificity matrix |
| 2 | **conLDSC-GES** | SC-to-GWAS | CELLECT-LDSC on the CELLEX GES matrix |
| 3 | **seismic-mBAT-combo** | GWAS-to-SC | `seismicGWAS` with mBAT-combo z-statistics substituted for MAGMA |
| 4 | **scDRS** | GWAS-to-SC | scDRS scored on an mBAT-combo-derived gene set |

MAGMA-GSEA and binary LDSC are benchmarked in the paper but are **not** CATCH components.
They remain in `modules/legacy/`, off by default — see `modules/legacy/README.md`.

## Pipeline structure

```
PREPARE_GENE_COORDS ──┬─> NORMALIZE_H5AD ──┬─> METRICS ──> CONLDSC (x2: cepo, ges) ─┐
                      │                     │                                        │
                      │                     ├─────────────────> SEISMIC ─────────────┤
                      │                     │                      ^                 ├─> COMBINE_CAUCHY
                      │                     └─────────────────> SCDRS ───────────────┤
                      │                                            ^                 │
                      └────────────> MBAT (gcta64 --mBAT-combo) ───┘                 │
                                                                                     v
                                                        <gwas>_catch_combined.tsv
```

`MBAT` runs once and is shared by the seismic and scDRS branches.

| Path | Purpose |
|---|---|
| `subworkflows/metrics.nf` | Cepo and CELLEX/GES specificity matrices |
| `subworkflows/conldsc.nf` | Native Nextflow port of CELLECT-LDSC, run once per specificity matrix |
| `subworkflows/mbat.nf` | GCTA mBAT-combo, shared |
| `subworkflows/seismic.nf` | mBAT z-statistics → `seismicGWAS` |
| `subworkflows/scdrs.nf` | mBAT gene set → scDRS |
| `modules/legacy/` | MAGMA-GSEA and the superseded LDSC implementation |

## Method details

These follow the manuscript's Methods; deviations are called out explicitly.

**Preprocessing.** `NORMALIZE_H5AD` reverses any prior log normalization, applies
UMI-based normalization, and emits **log2(TPM + 1)**. Cell types with **fewer than 20
cells** are dropped. Genes are restricted to **protein-coding autosomal** genes
(~19,430, matching CELLECT's `gene_coordinates.GRCh37.ensembl_v91.txt`). Cell type
labels are sanitized **once, here**, to CELLECT's `[A-Za-z0-9_-]` alphabet with no `__`,
so all four branches share identical labels and the Cauchy join cannot silently drop
cell types. The original→sanitized mapping is published alongside the h5ad.

**conLDSC.** A rule-for-rule port of CELLECT-LDSC (`perslab/CELLECT`,
`cellect-ldsc.snakefile`), not a reimplementation:

- ±100 kb around the gene body, clipped to chromosome size
- `bedops --partition` + `bedmap --echo-map-id-uniq` for disjoint overlap segments
- a SNP spanned by multiple genes gets the **maximum** ES value
- thin annotations, one combined file per chromosome, columns sorted alphabetically
- an **`all_genes_in_dataset` control annotation** alongside the 53-annotation baseline
- one `ldsc.py --h2-cts` job per GWAS across all cell types, via a `.ldcts` file
- the reported p-value is `Coefficient_P_value` = `norm.sf(coef/coef_se)`, the
  **one-sided coefficient z-score P value** the paper specifies (it states this
  outperforms heritability-enrichment p-values)

Output matches CELLECT's `prioritization.csv` schema — `gwas, specificity_id,
annotation, beta, beta_se, pvalue` — so it can be diffed column-for-column against a
reference CELLECT run.

Requires the **`pascaltimshel/ldsc` fork** pinned at
`d869cfd1e9fe1abc03b65c00b8a672bd530d0617`, not stock LDSC: its `--h2-cts` adds
per-cell-type error handling. Built by `environments/ldsc-timshel.def`.

**seismic-mBAT-combo.** Per the Methods, *"the disease genes z-statistics could be
replaced with mBAT-combo based z-statistics"*. `seismicGWAS::calc_specificity()`
computes its **own** specificity from the SingleCellExperiment — it does not consume
Cepo or GES, which is what keeps the four components independent. The z-statistic is
`qnorm(P_mBATcombo, lower.tail = FALSE)`: **one-sided**, matching MAGMA's `ZSTAT`
convention, because `get_ct_trait_associations` runs a one-sided test. Unlike the scDRS
branch there is **no top-N truncation** — seismic regresses across all overlapping genes,
and truncating would destroy the null part of the regression.

`seismicGWAS` has no CRAN/Bioconductor release and no git tags, so
`environments/seismic.def` pins it by commit.

**scDRS.** mBAT-combo p-values → top 1,000 genes → z-score weights → `scdrs munge-gs`
→ `compute-score` → `perform-downstream --group-analysis`. The per-cell-type p-value is
`assoc_mcp`. `--flag-raw-count` defaults to `False` because preprocessing emits
log2(TPM+1).

**Cauchy combination.** ACAT across the four component p-values **within each cell
type**. If the four inputs do not cover the same cell types the step **fails loudly**
with a set difference rather than inner-joining them away.

> **FDR is cross-trait.** The paper controls FDR at 5% "across all tissues/cell types
> **and traits** within each dataset". One pipeline run covers one trait, so
> `within_run_fdr` in the output is provisional. Run `bin/catch_fdr.R` over the
> `*_catch_combined.tsv` of all traits to reproduce the published thresholds.

## Requirements

- **Nextflow** >= 21.04.0
- **Singularity**
- A high-memory node: Cepo with `computePvalue = 100` and CELLEX both need ~256 GB

### Containers

Build these into `environments/` before the first run:

| Definition | Provides |
|---|---|
| `environments/ldsc-timshel.def` | `pascaltimshel/ldsc` @ d869cfd (`ldsc.py`, `munge_sumstats.py`) |
| `environments/cellect-py3.def` | pandas/pybedtools/bedtools/**bedops 2.4.37** for the CELLECT port |
| `environments/cellex.def` | CELLEX |
| `environments/seismic.def` | `seismicGWAS` pinned by commit |
| `py-r-cepo-scdrs.sif` | R stack (Cepo, zellkonverter, data.table) + scDRS |
| `gcta_v1.94.1.sif` | GCTA (`--mBAT-combo`) |

```bash
cd environments
for d in ldsc-timshel cellect-py3 cellex seismic; do
    singularity build --fakeroot $d.sif $d.def
done
```

### Reference data (hg19/GRCh37 throughout)

| Parameter | File |
|---|---|
| `ref_hg19_plink_prefix` | 1000G EUR Phase 3 PLINK, per chromosome |
| `ref_hg19_baseline` | Baseline annotations, **53 annotations**, thin-annot |
| `ref_hg19_weights` | `weights.hm3_noMHC.` |
| `ref_hg19_print_snps` | CELLECT `print_snps.txt` (1,217,311 rsIDs) |
| `ref_hg19_chr_sizes` | `GRCh37-chr-sizes.txt`, chr 1–22 |
| `ref_hg19_w_hm3_snplist` | `w_hm3.snplist`, for `munge_sumstats.py --merge-alleles` |
| `gene_matrix` | `geneMatrix.tsv.gz` |

Two easy mistakes:

- `--print-snps` must be CELLECT's `print_snps.txt`, **not** `hm3_no_MHC.list.txt`.
- The baseline must be the **53-annotation** model, not baselineLD v2.x (97 annotations).

No liftOver step is needed. CELLECT-LDSC and mBAT-combo are hg19 end to end,
`seismicGWAS` never reads coordinates (it joins on gene ID), and CELLEX is
coordinate-free. GWAS SNPs are reconciled to the reference by rsID.

### Input h5ad

- counts or log-normalized values in `adata.X` (declare which via `--h5ad_input_scale`)
- Ensembl gene IDs in `var_names`
- cell type labels in `adata.obs[<--cell_type_col>]`

## Installation

```bash
git clone https://github.com/PeterCAllen/CATCH-pipeline.git
cd CATCH-pipeline

curl -s https://get.nextflow.io | bash
mv nextflow ~/bin/

export SINGULARITY_CACHEDIR=$HOME/.singularity
mkdir -p $SINGULARITY_CACHEDIR
```

## Usage

```bash
nextflow run main.nf \
    --h5ad_input data/singlecell/human_atlas.h5ad \
    --gwas_sumstats data/gwas/AD_Jansen2019.txt.gz \
    --gwas_name AD_Jansen2019 \
    --genome_build hg19 \
    --gwas_sample_size 455258 \
    --gene_matrix data/reference/geneMatrix.tsv.gz \
    --cell_type_col cell.labels \
    --outdir results/human_atlas/AD_Jansen2019 \
    -profile pbs
```

`--gwas_name` must not contain `__` — that is CELLECT's `<specificity_id>__<annotation>`
separator.

`nextflow run main.nf --help` lists every parameter.

### Running components individually

Each component can be validated on its own before trusting the combination:

```bash
# conLDSC only
nextflow run main.nf --run_seismic false --run_scdrs false ...

# seismic and scDRS only (skips the expensive LDSC branch)
nextflow run main.nf --run_conldsc false ...

# legacy MAGMA, to validate seismic against its native MAGMA input
nextflow run main.nf --run_magma true --run_conldsc false --run_seismic false --run_scdrs false ...
```

`COMBINE_CAUCHY` runs only when all four components are enabled; otherwise the
per-component outputs are still produced and the pipeline warns.

### Cross-trait FDR

```bash
Rscript bin/catch_fdr.R results/catch_fdr.tsv results/*/*/combined/*_catch_combined.tsv
```

### Cluster configuration

Edit `conf/pbs.config` for your system:

```groovy
params {
    pbs_project = 'ab12'
    pbs_storage = 'scratch/ab12+gdata/ab12'
}
```

The queue is derived automatically from `task.memory`, so new processes should set a
`label` (`low_mem` / `medium_mem` / `high_mem`) or a `withName` block rather than
setting `queue` directly.

There is no SLURM profile.

## Output structure

```
<outdir>/
├── reference/                       gene coordinate files
├── preprocessed/                    normalized h5ad, cell type mapping and counts
├── metrics/
│   ├── cepo/                        cepo_norm.csv, cepo_pvalues.csv, cepo_s.csv
│   ├── cellex/                      ges.csv and the other CELLEX metrics
│   └── sanitized/                   CELLECT-ready matrices + annotation lists
├── conldsc/
│   ├── gwas/                        harmonized and munged sumstats
│   ├── cepo_norm/
│   │   ├── precomputation/          all_genes.*.csv, *.ldcts.txt
│   │   ├── control/                 all-genes control LD scores
│   │   ├── prioritization/          *.cell_type_results.txt
│   │   └── results/prioritization.csv
│   └── ges/                         same layout
├── mbat/                            per-chromosome and combined mBAT-combo results
├── seismic/                         z-statistics, specificity, *_seismic.tsv
├── scdrs/                           gene sets, cell scores, group results
└── combined/
    ├── <gwas>_catch_combined.tsv    CATCH p-values per cell type
    └── <gwas>_catch_plot.png
```

### Key output files

| File | Contents |
|---|---|
| `conldsc/<id>/results/prioritization.csv` | `gwas, specificity_id, annotation, beta, beta_se, pvalue` — CELLECT schema |
| `seismic/<gwas>_seismic.tsv` | `cell_type, pvalue, FDR` |
| `scdrs/results/*.scdrs_group.*` | per-cell-type `assoc_mcp`, `assoc_mcz`, `n_cell` |
| `combined/<gwas>_catch_combined.tsv` | per cell type: the four component p-values, `CATCH_P`, `within_run_fdr` |

## Validating against a reference run

Check in this order — a mismatch at any stage is far easier to debug than a wrong final
p-value:

1. **Preprocessing** — cell type list after the ≥20-cell filter; that everything
   downstream inherits it.
2. **Metrics** — `ges.csv` should be numerically *identical* to a CELLEX run on the same
   input. Cepo uses 100 permutations (seeded at 602), so expect small differences there.
3. **conLDSC** — compare in order: `genes_plus_100kb.<chr>.bed` row counts →
   `COMBINED_ANNOT.<chr>.annot.gz` column sums (these become `.l2.M`) →
   `.l2.ldscore.gz` correlation → `prioritization.csv`. Residual differences usually come
   from the baseline version, CELLECT's 1-based→0-based BED start offset, or pandas float
   formatting in the annot file.
4. **seismic** — no published reference exists for this pairing. Check that
   `qnorm(P_mBATcombo, lower.tail = FALSE)` is ~N(0,1) genome-wide, that gene overlap
   exceeds 80%, and cross-check against seismic run on MAGMA output via the legacy module.
5. **scDRS** — unchanged from previous releases apart from `--cov-file`,
   `--flag-return-ctrl-norm-score` and `--flag-raw-count False`.
6. **Cauchy** — the merged table must have the same cell type count as every input.

## Troubleshooting

**Out of memory.** Cepo and CELLEX both need ~256 GB. Adjust `conf/pbs.config` and use
`-resume`.

**`ldsc.py --h2-cts` fails or drops cell types.** LDSC hardcodes `_N_CHR = 22`, so all 22
chromosomes must have succeeded. Dropped cell types are logged as `*CTS ERROR*` and are
usually zero-variance annotations; `SANITIZE_SPECIFICITY` removes all-zero columns before
they reach LDSC.

**Cauchy step fails with a cell type set difference.** A branch dropped cell types.
Compare each component's output against `preprocessed/*.celltype_counts.tsv` to find
which one and why — do not work around it by loosening the join.

**`.gs` trait name is wrong.** `scdrs munge-gs` names the trait after the z-score column,
so `MBAT_TO_TSV` names that column after `--gwas_name`. `MUNGE_SCDRS_GENESET` asserts this.

**Missing reference files.** Only hg19 references are needed. Verify per-chromosome PLINK
files for chr 1–22 and that `print_snps.txt` and `GRCh37-chr-sizes.txt` are present.

```bash
nextflow run main.nf --help
```

## Citation

If you use this pipeline, please cite the CATCH paper along with the component methods:

- **CATCH**: Li A, Allen P, et al. *Benchmarking methods integrating GWAS and single-cell
  transcriptomic data for mapping trait-cell type associations.*
- **LDSC**: Bulik-Sullivan et al., Nature Genetics (2015)
- **CELLECT**: Timshel et al., eLife (2020)
- **Cepo**: Kim et al., Nature Computational Science (2021)
- **CELLEX**: Timshel et al., eLife (2020)
- **seismic**: Lai et al., Nature Communications (2025)
- **mBAT-combo**: Li et al., American Journal of Human Genetics (2023)
- **scDRS**: Zhang et al., Nature Genetics (2022)

## License

MIT License.

## Authors

- Ang Li
- Jian Zeng
- Peter C Allen

## Contact

https://github.com/PeterCAllen/CATCH-pipeline/issues
