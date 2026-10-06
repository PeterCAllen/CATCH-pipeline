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
- hg19/GRCh37 reference data under `data/reference/`. Fetch it all in one step, from a
  node with internet access (PBS compute nodes typically don't have one):

  ```bash
  bin/setup_reference_data.sh
  ```

  This pulls everything (1000G EUR Phase 3 PLINK panels, the baseline v1.1 LD score
  model, CELLECT's support files, the conLDSC gene-coordinate reference, `geneMatrix.tsv.gz`,
  and MAGMA's `g1000_eur` panel for mBAT-combo) from one source -- the
  ["CATCH Pipeline Reference" Zenodo record](https://doi.org/10.5281/zenodo.23174917)
  (CC-BY-4.0; see `MANIFEST.md` in that record for per-file provenance and citations --
  none of it is original work, it's a repackaging of public reference data). Verifies
  checksums and is safe to re-run (skips anything already present). Takes an optional
  target directory argument if you don't want the default `data/reference/`.

  `data/` is gitignored -- keep your copy on non-purged storage (e.g. `/g/data` on NCI
  systems), not `/scratch`.

Build the containers first (requires root or `--fakeroot`/subuid access -- on clusters
without that, build on your own machine with Docker/root and copy the resulting `.sif`
files onto a non-purged path such as `/g/data`, not `/scratch`):

```bash
cd environments
for d in ldsc-py3 cellect-py3 cellex seismic py-r-cepo-scdrs gcta_v1.94.1 scdrs-downstream; do
    singularity build --fakeroot $d.sif $d.def
done
```

`py-r-cepo-scdrs.def`, `gcta_v1.94.1.def`, and `scdrs-downstream.def` were authored by
inspecting how the pipeline calls into each image (see the `%labels` in each file);
`gcta_v1.94.1.def` pins an exact upstream version, and both `scdrs-downstream.def`
and the Python side of `py-r-cepo-scdrs.def` pin exact versions matching Ang Li's
environment (Scanpy 1.9.3, AnnData 0.8.0, NumPy 1.26.4, pandas 2.3.3, SciPy 1.14.1,
h5py 3.15.1, scDRS 1.0.4) for the Ang-vs-Peter CATCH benchmark -- see each file's
`%labels` for the source. `py-r-cepo-scdrs.def`'s R package versions are still
best-effort and were not verified against a prior working image -- check results
against any known-good reference before trusting them for publication numbers.

## Usage

```bash
nextflow run main.nf \
    --h5ad_input data/sc.h5ad \
    --gwas_sumstats data/gwas.txt.gz \
    --gwas_name AD_Jansen2019 \
    --genome_build hg19 \
    --gwas_sample_size 455258 \
    --gene_matrix data/reference/geneMatrix.tsv.gz \
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
| `bin/setup_reference_data.sh` | fetches all hg19/GRCh37 reference data (see Requirements) |
| `tests/` | tests for the CELLECT ports |
| `validation/ang_benchmark/` | cross-validation against an independent reimplementation of this analysis (not part of the pipeline itself) |

## Citation

Cite the CATCH paper along with the component methods: LDSC (Bulik-Sullivan 2015),
CELLECT and CELLEX (Timshel 2020), Cepo (Kim 2021), seismic (Lai 2025),
mBAT-combo (Li 2023), scDRS (Zhang 2022).

## Authors

Ang Li, Jian Zeng, Peter C Allen

Issues: https://github.com/PeterCAllen/CATCH-pipeline/issues
