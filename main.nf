#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Import modules and subworkflows
include { PREPARE_GENE_COORDS } from './modules/prepare_gene_coords'
include { RUN_CEPO            } from './modules/run_cepo'
include { MAGMA               } from './subworkflows/magma'
include { LDSC                } from './subworkflows/ldsc'
include { SCDRS               } from './subworkflows/scdrs'
include { COMBINE_CAUCHY      } from './modules/combine_cauchy'

// Validate required parameters
def validateParams() {
    def errors = []
    
    if (!params.h5ad_input) errors.add("--h5ad_input is required")
    if (!params.gwas_sumstats) errors.add("--gwas_sumstats is required")
    if (!params.genome_build) errors.add("--genome_build is required (hg19/GRCh37 or hg38/GRCh38)")
    if (!params.gwas_sample_size) errors.add("--gwas_sample_size is required")
    if (!params.gene_matrix) errors.add("--gene_matrix is required")
    
    if (params.run_magma) {
        if (!params.ref_1000g_prefix) errors.add("--ref_1000g_prefix is required for MAGMA")
    }
    
    if (errors.size() > 0) {
        log.error "Parameter validation failed:"
        errors.each { log.error "  - ${it}" }
        log.error "\nRun with --help for usage information"
        System.exit(1)
    }
}

workflow {
    // Validate parameters
    validateParams()
    
    // Normalize genome build
    def genome_build = params.genome_build.toLowerCase()
    if (genome_build == 'grch37') genome_build = 'hg19'
    if (genome_build == 'grch38') genome_build = 'hg38'
    
    // --- 0. Log parameters and start pipeline ---
    log.info """
    =====================================================
     CATCHY Pipeline
    =====================================================
     Input Parameters:
     ----------------
     h5ad_input       : ${params.h5ad_input}
     gwas_sumstats    : ${params.gwas_sumstats}
     genome_build     : ${genome_build}
     gwas_sample_size : ${params.gwas_sample_size}
     cell_type_col    : ${params.cell_type_col}
     
     Analysis Options:
     -----------------
     run_ldsc         : ${params.run_ldsc}
     run_magma        : ${params.run_magma}
     run_scdrs        : ${params.run_scdrs}
     ldsc_window_kb   : ${params.ldsc_window_kb}
     magma_window_kb  : ${params.magma_window_kb}
     
     Output:
     -------
     outdir           : ${params.outdir}
    =====================================================
    """

    // --- 1. Input Channel Creation ---
    ch_h5ad          = Channel.fromPath(params.h5ad_input, checkIfExists: true)
    ch_gwas          = Channel.fromPath(params.gwas_sumstats, checkIfExists: true)
    ch_gene_matrix   = Channel.fromPath(params.gene_matrix, checkIfExists: true)
    ch_genome_build  = Channel.value(genome_build)
    ch_magma_bin     = Channel.fromPath(params.magma_bin, checkIfExists: true)
    
    // Determine reference prefix based on genome build
    def ref_prefix = (genome_build in ['hg38', 'grch38']) ?
        params.ref_hg38_plink_prefix : params.ref_hg19_plink_prefix
    
    // Create reference file channel for MAGMA
    // Collect all per-chromosome reference files (.bed, .bim, .fam for chr 1-22)
    // PLUS the combined files (without chromosome number)
    ch_magma_ref = Channel.from(1..22)
        .flatMap { chr -> 
            [
                file("${ref_prefix}.${chr}.bed", checkIfExists: true),
                file("${ref_prefix}.${chr}.bim", checkIfExists: true),
                file("${ref_prefix}.${chr}.fam", checkIfExists: true)
            ]
        }
        .mix(
            Channel.fromPath([
                "${ref_prefix}.bed",
                "${ref_prefix}.bim",
                "${ref_prefix}.fam"
            ], checkIfExists: true)
        )
        .collect()
    
    // --- 2. Prepare Gene Coordinates (genome-build specific) ---
    PREPARE_GENE_COORDS(
        ch_gene_matrix,
        ch_genome_build
    )

    // --- 3. Run CEPO Analysis ---
    RUN_CEPO(
        ch_h5ad,
        PREPARE_GENE_COORDS.out.cepo_coords
    )

    // --- 4. Run MAGMA Subworkflow (if enabled) ---
    if (params.run_magma) {
        MAGMA(
            ch_gwas,
            RUN_CEPO.out.stats_tsv,
            ch_magma_bin,
            PREPARE_GENE_COORDS.out.magma_loc,
            ch_magma_ref,
            ch_genome_build
        )
    }

    // --- 5. Run LDSC Subworkflow (if enabled) ---
    if (params.run_ldsc) {
        LDSC(
            ch_gwas,
            RUN_CEPO.out.stats_tsv,
            PREPARE_GENE_COORDS.out.ldsc_loc,
            ch_genome_build
        )
    }

    // --- 6. Run scDRS Subworkflow (if enabled) ---
    if (params.run_scdrs) {
        SCDRS(
            ch_gwas,
            ch_h5ad,
            PREPARE_GENE_COORDS.out.cepo_coords,
            ch_genome_build
        )
    }

    // --- 7. Combine P-values from All Methods (if all three are enabled) ---
    if (params.run_ldsc && params.run_magma && params.run_scdrs) {
        // Extract h5ad basename for dataset name
        def h5ad_name = file(params.h5ad_input).baseName.replaceAll(/\.h5ad$/, '')
        
        COMBINE_CAUCHY(
            LDSC.out.annotation_comparison,
            MAGMA.out.gsa_results.first(),
            SCDRS.out.group_results.flatten().filter { it.name.contains('.scdrs_group') }.first(),
            h5ad_name,
            params.gwas_name
        )
    }
}

workflow.onComplete {
    def combined_msg = (params.run_ldsc && params.run_magma && params.run_scdrs) ? 
        "\n     Combined results: ${params.outdir}/combined/${params.gwas_name}_cauchy_combined.tsv" : ""
    
    log.info """
    =====================================================
     Pipeline completed at: ${workflow.complete}
     Execution status: ${workflow.success ? 'SUCCESS' : 'FAILED'}
     Duration: ${workflow.duration}
     Results: ${params.outdir}${combined_msg}
    =====================================================
    """.stripIndent()
}
