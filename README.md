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

Prebuilt images of all seven containers, built from the recipes below, can be downloaded from https://doi.org/10.5281/zenodo.23905323; put the `.sif` files in
`environments/`. To build them yourself instead (requires root or `--fakeroot`/subuid access -- on clusters
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
`%labels` for the source. `py-r-cepo-scdrs.def`'s R packages are not pinned one by one;
they are fixed by the Bioconductor 3.18 release and a dated CRAN snapshot (see table
below). To record the exact versions resolved in your build, run
`singularity exec environments/py-r-cepo-scdrs.sif Rscript -e 'sessionInfo()'`.

### Software versions

All dependencies run inside the Singularity containers built from `environments/*.def`;
the host only needs Linux (x86_64), Nextflow >= 21.04.0 and Singularity/Apptainer.

| Container | Base image | Pinned versions |
|---|---|---|
| `ldsc-py3` | condaforge/mambaforge:23.3.1-1 | Python 3.9; LDSC (CBIIT/ldsc @ `6c67395`) |
| `cellect-py3` | condaforge/mambaforge:23.3.1-1 | Python 3.6.9, pandas 0.24.2, NumPy 1.16.3, pybedtools 0.8.0, bedtools 2.27.1, BEDOPS 2.4.37 |
| `cellex` | condaforge/mambaforge:23.3.1-1 | Python 3.9, CELLEX 1.2.2, NumPy < 1.24 (Scanpy unpinned; 1.10.3 in our build) |
| `py-r-cepo-scdrs` | bioconductor/bioconductor_docker:RELEASE_3_18 (R 4.3.3) | Cepo 1.8.0 and other R packages from Bioconductor 3.18 / CRAN snapshot 2024-06-03; system Python 3.10 with Scanpy 1.9.3, AnnData 0.8.0, NumPy 1.26.4, pandas 2.3.3, SciPy 1.14.1, h5py 3.15.1, scDRS 1.0.4 (@ `644a5c7`) |
| `scdrs-downstream` | condaforge/mambaforge:23.3.1-1 | Python 3.11, NumPy 1.26.4, pandas 2.3.3, SciPy 1.14.1, AnnData 0.8.0, Scanpy 1.9.3, h5py 3.15.1, scDRS 1.0.4 (@ `644a5c7`) |
| `seismic` | bioconductor/bioconductor_docker:RELEASE_3_18 (R 4.3.3) | seismicGWAS (ylaboratory/seismic @ `1d29b5d`); CRAN snapshot 2024-03-01 |
| `gcta_v1.94.1` | debian:12-slim | GCTA 1.94.1 (static Linux binary) |

### Tested on

NCI Gadi (Rocky Linux 8.10, kernel 4.18.0, x86_64), PBS Pro, Nextflow 24.04.4,
SingularityCE 3.11.3 (containers built with SingularityCE 3.11.4). The CI workflow (`.github/workflows/tests.yml`) runs only the CELLECT port tests on
Python 3.11 / ubuntu-latest.

### Install time

Building the seven containers took about 40 min on an Apple Silicon MacBook (Colima VM,
4 CPUs, 8 GB RAM, x86_64 emulated through Rosetta) and gives about 5 GB of `.sif` files.
`bin/setup_reference_data.sh` downloads 1.15 GB (5.5 GB once extracted); this took 38 min
on an NCI Gadi login node, limited by Zenodo's download speed. A full analysis needs a
~256 GB RAM node (Cepo/CELLEX), so it is not practical on a desktop computer.

## Demo

The `test` profile runs CATCH on a small demo set: 2,400 cells (8 cell types x 300 cells,
sampled evenly across the 9 donors) from the human fetal bone marrow 10x data of
Jardine et al. 2021, against the Jansen et al. 2019 Alzheimer's disease GWAS. It needs
the containers and reference data from Requirements.

1. Get the demo inputs into `data/demo/`:
   - `fbm_demo.h5ad`: download from https://doi.org/10.5281/zenodo.23905323, or rebuild it from
     `fig1b_fbm_scaled_gex_updated_dr_20210104.h5ad` (from https://fbm.cellatlas.io/):

     ```bash
     python bin/subset_h5ad.py \
         --input_h5ad fig1b_fbm_scaled_gex_updated_dr_20210104.h5ad \
         --output_h5ad data/demo/fbm_demo.h5ad \
         --cell_type_col cell.labels --donor_col orig.ident --cells_per_type 300 \
         --gene_id_col gene_ids-1 \
         --cell_types "CD14 monocyte" "stromal macrophage" "DC2" "neutrophil" \
                      "naive B cell" "mid erythroid" "CD4 T cell" "sinusoidal EC"
     ```
   - `AD_sumstats_Jansenetal_2019sept.txt.gz` (686 MB, md5
     `12439e06b6bc1d8310d62b626df0b024`) from GWAS Catalog GCST007320:
     https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007320/

2. Run:

   ```bash
   nextflow run main.nf -profile pbs,test --pbs_project <project> --pbs_storage <storage>
   ```

3. Check `results/demo/combined/AD_Jansen2019_catch_combined.tsv`: one row per cell type
   (8 rows) with columns `group`, `dataset`, `trait`, `conLDSC_Cepo_pval`, `scDRS_pval`,
   `CATCH_P` and `within_run_fdr`, sorted by `CATCH_P`. Compare it with
   `tests/demo_expected.tsv`. The myeloid cell types should rank top (stromal
   macrophage, CATCH_P ≈ 0.002; CD14 monocyte ≈ 0.009; neutrophil ≈ 0.015), with
   lymphoid, erythroid and endothelial cells null:

   ```
   group               conLDSC_Cepo_pval  scDRS_pval  CATCH_P
   stromal_macrophage  0.0561             0.000999    0.00196
   CD14_monocyte       0.00464            0.156       0.00904
   neutrophil          0.00755            0.266       0.0148
   DC2                 0.0432             0.847       0.113
   mid_erythroid       0.208              0.500       0.316
   naive_B_cell        0.214              0.597       0.360
   sinusoidal_EC       0.262              0.661       0.441
   CD4_T_cell          0.643              0.993       0.986
   ```

Expected run time: about 35 min of wall time on NCI Gadi (`-profile pbs,test`, from
scratch, short queue waits), 172 tasks using about 37 allocated core-hours in total. The
longest steps are the 22 per-chromosome mBAT-combo jobs (up to ~20 min each, run in
parallel) and the conLDSC LD score jobs (up to ~5 min each); peak memory of any single
step is 7.2 GB (mBAT-combo). Fewer cells only shortens the expression steps
(normalisation, Cepo, scDRS); mBAT-combo and the conLDSC LD score steps still run
genome-wide, so on a single 8-core machine without a scheduler expect roughly 5 hours
(estimated from the core-hours above, not measured).
Repeat runs reproduce `tests/demo_expected.tsv` exactly. The demo runs scDRS without a
covariate file.

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

The input h5ad must use **Ensembl gene IDs** (e.g. `ENSG00000141510`) as `var_names`:
genes are matched to the gene coordinate files by ID, and the run stops at
`NORMALIZE_H5AD` with "No h5ad var_names matched the gene coordinate file" if they are
gene symbols. If your object has symbols as `var_names` and the IDs in a `var` column,
swap them first, e.g. in Python:

```python
adata.var["gene_symbol"] = adata.var_names
adata.var_names = adata.var["gene_ids"].astype(str)   # your Ensembl ID column
adata.var_names_make_unique()
```

(`bin/subset_h5ad.py --gene_id_col <column>` does the same while subsetting.) It also
needs cell type labels in the `obs` column given by `--cell_type_col`.

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

By default the mBAT-combo gene list and the hg19 conLDSC gene coordinate file cover all
genes. `--gene_coords_protein_coding_only true` restricts both to protein-coding genes
(about a third of the mBAT-combo gene tests). This is separate from
`--protein_coding_only`, which filters the expression matrix (and, for hg38 runs, the
conLDSC coordinates as well). scDRS and seismic only use
genes present in the (protein-coding) h5ad, so their results should not change. conLDSC
results can shift slightly: CELLECT gives an ES of 0 to coordinate-file genes that are
absent from the specificity matrix, so dropping non-coding genes changes some
annotation values.

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

## Citation

Cite the CATCH paper along with the component methods: LDSC (Bulik-Sullivan 2015),
CELLECT and CELLEX (Timshel 2020), Cepo (Kim 2021), seismic (Lai 2025),
mBAT-combo (Li 2023), scDRS (Zhang 2022).

## Authors

Ang Li, Jian Zeng, Peter C Allen

Issues: https://github.com/PeterCAllen/CATCH-pipeline/issues
