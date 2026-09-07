#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//
// CATCH — Cell-type Associations with Traits via CaucHy combination of strategies
// (Li et al.). Cauchy combination of two components:
//   conLDSC-Cepo, scDRS (mBAT-combo based)
//
// conLDSC-GES and seismic-mBAT-combo remain available as optional, off-by-default
// validation branches (see --run_cellex / --run_seismic) but are not part of CATCH.
//

include { asBool } from './lib/util'

include { PREPARE_GENE_COORDS } from './modules/prepare_gene_coords'
include { NORMALIZE_H5AD      } from './modules/normalize_h5ad'
include { CONVERT_GWAS_FOR_LDSC } from './modules/munge_gwas_for_ldsc'
include { MUNGE_SUMSTATS      } from './modules/munge_sumstats'
include { COMBINE_CAUCHY      } from './modules/combine_cauchy'

include { METRICS  } from './subworkflows/metrics'
include { CONLDSC  } from './subworkflows/conldsc'
include { MBAT     } from './subworkflows/mbat'
include { SEISMIC  } from './subworkflows/seismic'
include { SCDRS    } from './subworkflows/scdrs'

include { MAGMA as LEGACY_MAGMA } from './subworkflows/legacy/magma'

def helpMessage() {
    log.info """
    =====================================================
     ${workflow.manifest.name} v${workflow.manifest.version}
    =====================================================
    Usage:
      nextflow run main.nf --h5ad_input <file> --gwas_sumstats <file> --genome_build <hg19|hg38> --gwas_sample_size <N> [options]

    CATCH combines two components via Cauchy combination (Li et al.):
      conLDSC-Cepo, scDRS (mBAT-combo based)

    Required Arguments:
      --h5ad_input PATH          Path to input h5ad file
      --gwas_name STR            Short trait identifier (must not contain '__')
      --genome_build STR         Genome build (hg19/GRCh37 or hg38/GRCh38)
      --gwas_sample_size INT     Total sample size for GWAS study

    GWAS summary statistics (at least one required):
      --gwas_sumstats PATH       Raw GWAS sumstats (gzipped); reformatted internally
                                 for each branch that needs it
      --gwas_cojo PATH           Pre-formatted GCTA-COJO .ma (SNP A1 A2 freq b se p N);
                                 bypasses reformatting for the seismic/scDRS branch
      --gwas_sumstats_munged PATH  Pre-munged LDSC .sumstats.gz; bypasses munging for
                                 the conLDSC branch
      --gene_matrix PATH         Path to geneMatrix.tsv.gz file
      --cell_type_col STR        Column in h5ad.obs with cell types [${params.cell_type_col}]

    Preprocessing (Li et al.: log2(TPM+1), protein-coding, >=20 cells per cell type):
      --protein_coding_only BOOL    Restrict to protein-coding genes [${params.protein_coding_only}]
      --min_cells_per_celltype INT  Drop cell types below this [${params.min_cells_per_celltype}]
      --h5ad_input_scale STR        auto | raw_counts | lognorm [${params.h5ad_input_scale}]
      --h5ad_log_base STR           Base of any prior log normalization: e | 2 [${params.h5ad_log_base}]
      --h5ad_use_raw BOOL           Promote adata.raw to X first [${params.h5ad_use_raw}]

    Specificity metrics:
      --run_cepo BOOL               Compute Cepo [${params.run_cepo}]
      --run_cellex BOOL             Compute CELLEX (provides GES); optional validation
                                     branch, not part of CATCH [${params.run_cellex}]
      --cepo_compute_pvalue INT     Cepo permutations [${params.cepo_compute_pvalue}]
      --cepo_prefilter_pzero FLOAT  Cepo prefilter_pzero [${params.cepo_prefilter_pzero}]

    conLDSC (CELLECT-LDSC):
      --run_conldsc BOOL            Run conLDSC [${params.run_conldsc}]
      --conldsc_cepo_variant STR    cepo_norm | cepo_s [${params.conldsc_cepo_variant}]
      --conldsc_window_kb INT       Gene window in kb [${params.conldsc_window_kb}]
      --conldsc_ld_wind_cm FLOAT    LD window in cM [${params.conldsc_ld_wind_cm}]
      --conldsc_top_pct FLOAT       Restrict to top fraction of genes per cell type;
                                    null uses the full ES matrix as CELLECT does [${params.conldsc_top_pct}]
      --conldsc_keep_annots BOOL    Also write full SNP:ES annot files [${params.conldsc_keep_annots}]

    mBAT-combo (shared by seismic and scDRS):
      --mbat_window_kb INT          TSS/TES window in kb [${params.mbat_window_kb}]

    seismic (optional validation branch, not part of CATCH):
      --run_seismic BOOL            Run seismic-mBAT-combo [${params.run_seismic}]
      --seismic_assay STR           Assay for calc_specificity [${params.seismic_assay}]
      --seismic_lognorm BOOL        Run scater::logNormCounts first [${params.seismic_lognorm}]

    scDRS:
      --run_scdrs BOOL              Run scDRS [${params.run_scdrs}]
      --scdrs_top_genes INT         Top genes for the gene set [${params.scdrs_top_genes}]
      --scdrs_n_ctrl INT            Control gene sets [${params.scdrs_n_ctrl}]
      --scdrs_filter_data STR       Filter data ("True"/"False") [${params.scdrs_filter_data}]
      --scdrs_raw_count STR         Raw counts ("True"/"False") [${params.scdrs_raw_count}]
      --scdrs_cov_file PATH         Covariate file [${params.scdrs_cov_file}]

    Reference Data (hg19/GRCh37 throughout):
      --ref_hg19_plink_prefix PATH  1000G EUR Phase 3 PLINK prefix [${params.ref_hg19_plink_prefix}]
      --ref_hg19_baseline PATH      Baseline prefix, 53 annotations [${params.ref_hg19_baseline}]
      --ref_hg19_weights PATH       LD weights prefix (hm3_noMHC) [${params.ref_hg19_weights}]
      --ref_hg19_print_snps PATH    CELLECT print_snps.txt [${params.ref_hg19_print_snps}]
      --ref_hg19_chr_sizes PATH     GRCh37-chr-sizes.txt [${params.ref_hg19_chr_sizes}]
      --ref_hg19_w_hm3_snplist PATH w_hm3.snplist for munge_sumstats [${params.ref_hg19_w_hm3_snplist}]

    Legacy (not CATCH components; see modules/legacy/README.md):
      --run_magma BOOL              Run MAGMA-GSEA [${params.run_magma}]
      --magma_bin PATH              Path to MAGMA binary [${params.magma_bin}]

    Profiles:
      -profile standard          Run locally with Singularity
      -profile pbs               Submit jobs to PBS cluster (NCI Gadi)

    Example:
      nextflow run main.nf \\
        --h5ad_input data/sc.h5ad \\
        --gwas_sumstats data/gwas.txt.gz \\
        --gwas_name AD_Jansen2019 \\
        --genome_build hg19 \\
        --gwas_sample_size 455258 \\
        --gene_matrix data/geneMatrix.tsv.gz \\
        --cell_type_col cell.labels \\
        -profile pbs

    =====================================================
    """.stripIndent()
}

def validateParams(flags) {
    def errors = []

    if (!params.h5ad_input)      errors.add("--h5ad_input is required")
    if (!params.gwas_name)       errors.add("--gwas_name is required")
    if (!params.genome_build)    errors.add("--genome_build is required (hg19/GRCh37 or hg38/GRCh38)")
    if (!params.gwas_sample_size) errors.add("--gwas_sample_size is required")
    if (!params.gene_matrix)     errors.add("--gene_matrix is required")
    if (!params.cell_type_col)   errors.add("--cell_type_col is required")

    if (!params.gwas_sumstats && !params.gwas_cojo && !params.gwas_sumstats_munged) {
        errors.add("at least one of --gwas_sumstats, --gwas_cojo, --gwas_sumstats_munged is required")
    }

    if (params.gwas_name?.contains('__')) {
        errors.add("--gwas_name must not contain '__' (the CELLECT '<id>__<annotation>' separator)")
    }

    if (flags.conldsc && !flags.cepo && !flags.cellex) {
        errors.add("--run_conldsc requires at least one of --run_cepo / --run_cellex")
    }

    if (flags.conldsc && !params.gwas_sumstats && !params.gwas_sumstats_munged) {
        errors.add("--run_conldsc requires --gwas_sumstats or --gwas_sumstats_munged")
    }

    if ((flags.seismic || flags.scdrs) && !params.gwas_sumstats && !params.gwas_cojo) {
        errors.add("--run_seismic/--run_scdrs require --gwas_sumstats or --gwas_cojo")
    }

    if (flags.magma) {
        if (!params.magma_bin) errors.add("--magma_bin is required when --run_magma is set")
        if (!flags.cepo)       errors.add("--run_magma requires --run_cepo (gene sets are built from Cepo)")
        if (!params.gwas_sumstats) errors.add("--run_magma requires --gwas_sumstats (no pre-formatted bypass for MAGMA)")
    }

    if (errors.size() > 0) {
        log.error "Parameter validation failed:"
        errors.each { log.error "  - ${it}" }
        System.exit(1)
    }
}

workflow {
    if (params.help) {
        helpMessage()
        return
    }

    def flags = [
        cepo:    asBool(params.run_cepo),
        cellex:  asBool(params.run_cellex),
        conldsc: asBool(params.run_conldsc),
        seismic: asBool(params.run_seismic),
        scdrs:   asBool(params.run_scdrs),
        magma:   asBool(params.run_magma),
    ]

    validateParams(flags)

    def genome_build = params.genome_build.toLowerCase()
    if (genome_build == 'grch37') genome_build = 'hg19'
    if (genome_build == 'grch38') genome_build = 'hg38'

    def catch_ready = flags.conldsc && flags.cepo && flags.scdrs

    log.info """
    =====================================================
     CATCH Pipeline
    =====================================================
     Inputs
     ------
     h5ad_input       : ${params.h5ad_input}
     gwas_sumstats    : ${params.gwas_sumstats}
     gwas_name        : ${params.gwas_name}
     genome_build     : ${genome_build}
     gwas_sample_size : ${params.gwas_sample_size}
     cell_type_col    : ${params.cell_type_col}

     CATCH Components
     ----------------
     conLDSC-Cepo         : ${flags.conldsc && flags.cepo} (variant: ${params.conldsc_cepo_variant})
     scDRS                : ${flags.scdrs}
     CATCH ready          : ${catch_ready}

     Optional validation branches (not part of CATCH)
     --------------------------------------------------
     conLDSC-GES          : ${flags.conldsc && flags.cellex}
     seismic-mBAT-combo   : ${flags.seismic}
     MAGMA-GSEA           : ${flags.magma}

     Output
     ------
     outdir           : ${params.outdir}
    =====================================================
    """.stripIndent()

    ch_h5ad         = Channel.fromPath(params.h5ad_input,    checkIfExists: true)
    ch_gwas         = params.gwas_sumstats
        ? Channel.fromPath(params.gwas_sumstats, checkIfExists: true)
        : Channel.empty()
    ch_gene_matrix  = Channel.fromPath(params.gene_matrix,   checkIfExists: true)
    ch_genome_build = Channel.value(genome_build)

    // --- 1. Reference gene coordinates ---
    PREPARE_GENE_COORDS(ch_gene_matrix, ch_genome_build)

    // --- 2. Standardize the single-cell reference ---
    NORMALIZE_H5AD(ch_h5ad, PREPARE_GENE_COORDS.out.cepo_coords)
    ch_norm_h5ad = NORMALIZE_H5AD.out.h5ad

    ch_prioritization = Channel.empty()
    ch_seismic        = Channel.empty()
    ch_scdrs_group    = Channel.empty()
    ch_specificity    = Channel.empty()

    // --- 3. Specificity metrics (also needed by the legacy MAGMA gene sets) ---
    if (flags.conldsc || flags.magma) {
        METRICS(ch_norm_h5ad)
        ch_specificity = METRICS.out.specificity
    }

    // --- 3b. conLDSC-Cepo and conLDSC-GES ---
    if (flags.conldsc) {
        if (params.gwas_sumstats_munged) {
            // Pre-munged LDSC sumstats supplied directly; skip CONVERT_GWAS_FOR_LDSC + MUNGE_SUMSTATS.
            ch_gwas_munged = Channel.fromPath(params.gwas_sumstats_munged, checkIfExists: true)
        } else {
            ch_ref_bim = Channel.value(file(params.ref_hg19_bim_file, checkIfExists: true))
            CONVERT_GWAS_FOR_LDSC(ch_gwas, ch_genome_build, ch_ref_bim)

            MUNGE_SUMSTATS(
                CONVERT_GWAS_FOR_LDSC.out.formatted_gwas,
                Channel.value(file(params.ref_hg19_w_hm3_snplist, checkIfExists: true))
            )

            ch_gwas_munged = MUNGE_SUMSTATS.out.gwas_munged
        }

        CONLDSC(
            ch_specificity,
            METRICS.out.annotations,
            ch_gwas_munged,
            PREPARE_GENE_COORDS.out.cellect_loc
        )

        ch_prioritization = CONLDSC.out.prioritization
    }

    // --- 4. mBAT-combo, shared by seismic and scDRS ---
    if (flags.seismic || flags.scdrs) {
        MBAT(ch_gwas, PREPARE_GENE_COORDS.out.cepo_coords, ch_genome_build)

        if (flags.seismic) {
            SEISMIC(
                MBAT.out.mbat_combined,
                ch_norm_h5ad,
                PREPARE_GENE_COORDS.out.cepo_coords
            )
            ch_seismic = SEISMIC.out.results
        }

        if (flags.scdrs) {
            ch_scdrs_cov = params.scdrs_cov_file
                ? Channel.fromPath(params.scdrs_cov_file, checkIfExists: true)
                : Channel.value(file("NO_FILE"))

            SCDRS(MBAT.out.mbat_combined, ch_norm_h5ad, ch_scdrs_cov)
            ch_scdrs_group = SCDRS.out.group_results.first()
        }
    }

    // --- 5. Legacy MAGMA-GSEA (off by default; used to validate the seismic branch) ---
    if (flags.magma) {
        ch_magma_ref = Channel.from(1..22)
            .flatMap { chr ->
                ['bed', 'bim', 'fam'].collect { ext ->
                    file("${params.ref_hg19_plink_prefix}.${chr}.${ext}", checkIfExists: true)
                }
            }
            .mix(Channel.fromPath(
                ['bed', 'bim', 'fam'].collect { "${params.ref_hg19_plink_combined}.${it}" },
                checkIfExists: true))
            .collect()

        ch_magma_specificity = ch_specificity
            .filter { id, m -> id in ['cepo_norm', 'cepo_s'] }
            .map { id, m -> m }
            .first()

        LEGACY_MAGMA(
            ch_gwas,
            ch_magma_specificity,
            Channel.fromPath(params.magma_bin, checkIfExists: true),
            PREPARE_GENE_COORDS.out.magma_loc,
            ch_magma_ref,
            ch_genome_build
        )
    }

    // --- 6. CATCH: Cauchy combination of the two components ---
    if (catch_ready) {
        def h5ad_name = file(params.h5ad_input).simpleName

        ch_cepo_prior = ch_prioritization
            .filter { id, f -> id in ['cepo_norm', 'cepo_s'] }
            .map { id, f -> f }

        // Build the --method specs from the real staged filenames so they cannot drift.
        def CATCH_ORDER = ['conLDSC_Cepo', 'scDRS']

        ch_methods = ch_cepo_prior.map { f -> tuple('conLDSC_Cepo', "conLDSC_Cepo=${f.name}=annotation=pvalue", f) }
            .mix(ch_scdrs_group.map    { f -> tuple('scDRS',        "scDRS=${f.name}=group=assoc_mcp",          f) })
            .toList()
            .map { rows ->
                if (rows.size() != CATCH_ORDER.size()) {
                    error "CATCH expected ${CATCH_ORDER.size()} component result files, got ${rows.size()}: " +
                          rows.collect { it[0] }
                }
                def sorted = rows.sort { a, b -> CATCH_ORDER.indexOf(a[0]) <=> CATCH_ORDER.indexOf(b[0]) }
                tuple(sorted.collect { it[2] }, sorted.collect { it[1] })
            }

        COMBINE_CAUCHY(
            ch_methods,
            h5ad_name,
            params.gwas_name
        )
    }
    else {
        log.warn "CATCH combination skipped: --run_conldsc, --run_cepo, and --run_scdrs must all be enabled. " +
                 "Per-component results are still produced."
    }
}

