# CATCHY Pipeline

**C**ell-type **A**ssociation **T**esting with **C**EP**O** and **H**eritabilit**Y**

A Nextflow pipeline for integrating single-cell RNA-seq data with GWAS summary statistics using three complementary methods: LDSC (Linkage Disequilibrium Score Regression), MAGMA (Multi-marker Analysis of GenoMic Annotation), and scDRS (single-cell Disease Relevance Score).

## Overview

This pipeline performs cell-type-specific enrichment analysis by:

1. **Gene Coordinate Preparation**: Extracts hg19 gene coordinates from GeneMatrix (used by all analyses regardless of GWAS build)
2. **CEPO Analysis**: Identifies differentially expressed genes for each cell type using CEPO (Cell identity from Expression Profiles Optimization)
3. **MAGMA**: Gene set enrichment analysis using GWAS summary statistics
4. **LDSC**: Stratified LD score regression for heritability enrichment
5. **scDRS**: Disease relevance scoring at the single-cell level
6. **Cauchy Combination**: Combines p-values from all three methods

### Important: Reference Data Simplification

**All analyses (MAGMA, LDSC, scDRS) use hg19/GRCh37 reference data only.** If your GWAS is in hg38/GRCh38:
- **Gene coordinates** are automatically set to hg19 for all analyses (extracted from GeneMatrix)
- **SNP coordinates** in the GWAS file retain their original genome build (hg19 or hg38)
- All reference panels (1000G, LD scores) are hg19 only
- This approach works because:
  - Gene-SNP mapping uses the gene window around gene coordinates (hg19)
  - SNP matching with reference panels uses rsIDs, not coordinates
  - This eliminates the need for duplicate reference datasets in both genome builds

## Requirements

### Software Dependencies

- **Nextflow** (>=21.04.0)
- **Singularity** or **Docker** (for containerized execution)
- At least 256GB RAM for CEPO analysis (high-memory node recommended)

### Reference Data

**Important**: This pipeline always uses **hg19/GRCh37 reference data** for MAGMA, LDSC, and scDRS (mBAT). When your GWAS is in hg38/GRCh38, gene coordinates are automatically converted to hg19 using the GeneMatrix file, ensuring compatibility with reference data.

#### Required for all analyses:
- **Gene matrix file** (`geneMatrix.tsv.gz`): Gene coordinates for both hg19 and hg38 (used for automatic coordinate conversion)

#### For MAGMA and LDSC (hg19 only):
- **1000 Genomes Phase 3 European Panel (hg19) v1.1**:
  - PLINK binary files (`.bed`, `.bim`, `.fam`) per chromosome and combined
  - Example prefix: `/path/to/1000G.EUR.hg19`

#### For LDSC (hg19 only):
- **Baseline LD annotations**: `baselineLD.[1-22].annot.gz` and `.M` files
- **LD score weights**: `weights.hm3_noMHC.[1-22].l2.ldscore.gz`
- **HapMap3 SNP list**: `hm3_no_MHC.list.txt`

#### For MAGMA:
- **MAGMA binary**: Download from https://cncr.nl/research/magma/

### Data Format Requirements

#### Single-cell H5AD file:
- Must contain:
  - Raw or log-normalized counts in `adata.X`
  - Cell type annotations in `adata.obs` (column specified by `--cell_type_col`)
  - Gene names in `adata.var_names`

#### scDRS Covariate file (optional but recommended):
- **Format**: Tab-delimited (TSV) file
- **Required columns**:
  - `SAMPLE_ID`: Cell barcode matching those in the h5ad file
  - `N_GENE`: Number of genes detected per cell (technical covariate)
  - `SEX`: Biological sex (categorical: M/F or 0/1)
  - `AGE`: Age in years (numeric)
- **Purpose**: Controls for technical (sequencing depth) and biological (sex, age) confounders in scDRS disease association testing
- **Note**: This file should be prepared during single-cell pre-processing and quality control steps

#### GWAS Summary Statistics:
- Tab-delimited or space-delimited text file (gzipped)
- Must contain columns:
  - SNP ID or rsID
  - Chromosome
  - Base pair position
  - Effect allele (A1)
  - Other allele (A2)
  - P-value
  - Optional: Beta, SE, N (sample size)

## Installation

1. **Clone the repository**:
```bash
git clone https://github.com/PeterCAllen/CATCHY-pipeline.git
cd CATCHY-pipeline
```

2. **Install Nextflow** (if not already installed):
```bash
curl -s https://get.nextflow.io | bash
mv nextflow ~/bin/  # or add to your PATH
```

3. **Set up Singularity cache** (recommended):
```bash
export SINGULARITY_CACHEDIR=$HOME/.singularity
mkdir -p $SINGULARITY_CACHEDIR
```

## Configuration

### Required Parameters

Create a custom configuration file or specify parameters on the command line:

```bash
# Create a custom config file (e.g., my_config.config)
params {
    // ===== REQUIRED INPUT/OUTPUT =====
    h5ad_input       = '/path/to/your/singlecell_data.h5ad'
    gwas_sumstats    = '/path/to/your/gwas_summary_stats.txt.gz'
    gwas_name        = 'MyGWAS'  // Short name for output files
    outdir           = 'results/singlecell_data/MyGWAS'  // Recommended format
    
    // ===== GENOME BUILD (CRITICAL!) =====
    genome_build     = 'hg19'  // Options: hg19, GRCh37, hg38, GRCh38
    gwas_sample_size = 455258  // Total sample size of GWAS
    
    // ===== SINGLE-CELL PARAMETERS =====
    cell_type_col    = 'cell_type'  // Column name with cell type annotations
    
    // ===== GENE COORDINATES =====
    gene_matrix      = '/path/to/geneMatrix.tsv.gz'
    
    // ===== REFERENCE DATA: hg19/GRCh37 (ONLY) =====
    // Note: Pipeline uses hg19 references for all analyses
    // hg38 GWAS coordinates are automatically converted
    ref_hg19_plink_prefix = '/path/to/1000G.EUR.hg19'
    ref_hg19_baseline     = '/path/to/hg19/baseline/baselineLD.'
    ref_hg19_weights      = '/path/to/hg19/weights/weights.hm3_noMHC.'
    ref_hg19_hapmap3      = '/path/to/hg19/hm3_no_MHC.list.txt'
    
    // ===== MAGMA BINARY =====
    magma_bin        = '/path/to/magma'  // or just 'magma' if in PATH
    
    // ===== OPTIONAL: scDRS COVARIATES =====
    scdrs_cov_file   = '/path/to/covariates.tsv'  // Optional: TSV with SAMPLE_ID, N_GENE, SEX, AGE
    
    // ===== OPTIONAL: Analysis toggles =====
    run_magma        = true
    run_ldsc         = true
    run_scdrs        = true
}
```

### Output Directory Format

**Recommended format**: `results/<single-cell-object>/<gwas-name>`

Example:
```bash
--outdir results/PBMC_10k/Height_GWAS
--outdir results/brain_cortex/Alzheimers_2023
--outdir results/liver_tissue/Type2Diabetes
```

This structure helps organize results when analyzing:
- Multiple GWAS with the same single-cell dataset
- Multiple single-cell datasets with the same GWAS
- Many combinations of both

### Cluster Configuration

#### PBS Cluster

If running on a PBS cluster, edit `conf/pbs.config` to match your system:

```bash
params {
    pbs_project = 'ab12'                           // Your PBS project code
    pbs_storage = 'scratch/ab12+gdata/ab12'       // Your storage paths
}
```

You may also need to adjust:
- Queue names in the `process.queue` selector
- Module names if your cluster uses different environment modules
- Scratch directory paths

#### SLURM Cluster

If running on a SLURM cluster, edit `conf/slurm.config` to match your system:

```bash
params {
    slurm_account = 'def-username'                 // Your SLURM account/project
    slurm_partition = 'compute'                    // Default partition name
}
```

You may also need to adjust:
- Partition names in the `process.queue` selector (e.g., 'highmem', 'hugemem')
- Memory thresholds for automatic partition selection
- Add `clusterOptions` for any cluster-specific SLURM flags (e.g., QoS)
- Module loading if your cluster requires it

## Usage

### Basic Usage (Local Execution)

```bash
nextflow run main.nf \
  --h5ad_input data/pbmc_10k.h5ad \
  --gwas_sumstats data/height_gwas.txt.gz \
  --gwas_name Height \
  --genome_build hg19 \
  --gwas_sample_size 253288 \
  --gene_matrix data/reference/geneMatrix.tsv.gz \
  --cell_type_col cell_type \
  --outdir results/pbmc_10k/Height \
  --ref_hg19_plink_prefix data/reference/1000G.EUR.hg19 \
  --ref_hg19_baseline data/reference/hg19/baseline/baselineLD. \
  --ref_hg19_weights data/reference/hg19/weights/weights.hm3_noMHC. \
  --ref_hg19_hapmap3 data/reference/hg19/hm3_no_MHC.list.txt \
  --magma_bin resources/magma \
  --scdrs_cov_file data/pbmc_10k_covariates.tsv \
  -profile standard
```

### Cluster Execution

#### PBS Cluster

```bash
nextflow run main.nf \
  -c my_config.config \
  -profile pbs \
  -resume
```

#### SLURM Cluster

```bash
nextflow run main.nf \
  -c my_config.config \
  -profile slurm \
  -resume
```

### Using a Custom Config File

```bash
# PBS cluster
nextflow run main.nf -c my_config.config -profile pbs

# SLURM cluster
nextflow run main.nf -c my_config.config -profile slurm
```

## Analysis Options

### Toggle Individual Analyses

```bash
--run_magma true   # Enable/disable MAGMA (default: true)
--run_ldsc true    # Enable/disable LDSC (default: true)
--run_scdrs true   # Enable/disable scDRS (default: true)
```

### Window Sizes

```bash
--magma_window_kb 10    # MAGMA gene window (default: 10kb)
--ldsc_window_kb 100    # LDSC gene window (default: 100kb)
--scdrs_mbat_window_kb 50  # scDRS mBAT window (default: 50kb)
```

### scDRS-Specific Options

```bash
--scdrs_top_genes 1000      # Top genes for gene set (default: 1000)
--scdrs_n_ctrl 1000         # Control gene sets (default: 1000)
--scdrs_filter_data "False" # Filter cells/genes: "True" or "False"
--scdrs_raw_count "True"    # Expect raw counts: "True" or "False"
--scdrs_cov_file covariates.tsv  # Optional: Covariate file (TSV format)
```

**Important**: 
- `scdrs_filter_data` and `scdrs_raw_count` must be strings ("True" or "False"), not boolean values.
- The covariate file is optional but recommended for controlling technical and biological confounders
- Covariate file must contain columns: `SAMPLE_ID`, `N_GENE`, `SEX`, `AGE`

## Output Structure

```
results/<single-cell-object>/<gwas-name>/
├── cepo/
│   ├── <dataset>_cepo_stats.tsv           # CEPO statistics per cell type
│   └── <dataset>_cepo_stats.rds           # CEPO results (R object)
├── magma/
│   ├── <celltype>/
│   │   ├── magma_geneset.txt              # Gene set definition
│   │   ├── magma.genes.out                # Gene-level results
│   │   └── magma.gsa.out                  # Gene set enrichment results
│   └── magma_gsa_results.tsv              # Combined results across cell types
├── ldsc/
│   ├── <celltype>/
│   │   ├── ldsc_*.results                 # LDSC output files
│   │   └── ldsc_*.log                     # LDSC log files
│   └── annotation_comparison.tsv          # Combined LDSC results
├── scdrs/
│   ├── <dataset>.scdrs_group.<gwas>.tsv   # Cell type group results
│   ├── <dataset>.scdrs_ct.<gwas>.tsv      # Cell type results
│   └── <dataset>.full_score.<gwas>.gz     # Full per-cell scores
└── combined/
    └── <gwas>_cauchy_combined.tsv         # Combined p-values from all methods
```

## Key Output Files

### CEPO Results
- `*_cepo_stats.tsv`: Differentially expressed genes per cell type with statistics

### MAGMA Results
- `magma_gsa_results.tsv`: Cell type enrichment p-values and betas

### LDSC Results
- `annotation_comparison.tsv`: Heritability enrichment per cell type

### scDRS Results
- `*.scdrs_group.*.tsv`: Disease association scores per cell type
- `*.scdrs_ct.*.tsv`: Detailed cell type statistics

### Combined Results
- `*_cauchy_combined.tsv`: Integrated p-values using Cauchy combination test

## Troubleshooting

### Common Issues

1. **Out of memory errors**:
   - CEPO analysis requires high memory (256GB recommended)
   - Adjust `conf/pbs.config` resource allocations
   - Use `-resume` to restart from checkpoint

2. **Missing reference files**:
   - Verify all hg19 reference paths exist
   - Only hg19 references are required (hg38 coordinates are automatically converted)
   - Ensure per-chromosome PLINK files are present (chr 1-22)
   - Verify combined reference files exist (without chromosome suffix)

3. **Singularity/container issues**:
   - Set `SINGULARITY_CACHEDIR` to writable location
   - Check Singularity is available: `singularity --version`
   - Try rebuilding container cache

4. **Pipeline fails at specific step**:
   - Check `.nextflow.log` for detailed error messages
   - Examine `work/` directory for process-specific logs
   - Use `-resume` to continue from last successful step

### Getting Help

```bash
nextflow run main.nf --help
```

## Citation

If you use this pipeline, please cite:

- **LDSC**: Finucane et al., Nature Genetics (2015)
- **MAGMA**: de Leeuw et al., PLOS Computational Biology (2015)
- **scDRS**: Zhang et al., Nature Genetics (2022)
- **CEPO**: Wang et al., Nature Computational Science (2021)

## License

This pipeline is distributed under the MIT License.

## Authors

- Ang Li
- Zhen (Jennifer) Gao
- Jian Zeng
- Peter C Allen

## Contact

For questions or issues, please open an issue on the GitHub repository:
https://github.com/PeterCAllen/CATCHY-pipeline/issues
