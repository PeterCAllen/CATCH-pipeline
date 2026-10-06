// subworkflows/magma.nf

include { CREATE_MAGMA_GENESET } from '../../modules/legacy/create_magma_geneset'
include { RUN_MAGMA } from '../../modules/legacy/run_magma'

workflow MAGMA {
    take:
    gwas_ch           // path: GWAS summary statistics
    cepo_stats_ch     // path: CEPO statistics TSV
    magma_bin_ch      // path: MAGMA binary
    gene_loc_ch       // path: MAGMA gene location file
    ref_files_ch      // path: reference files (chr 1-22 .bed, .bim, .fam)
    genome_build_ch   // val: genome build
    
    main:
    // Step 1: Create MAGMA gene sets from CEPO results
    CREATE_MAGMA_GENESET(
        cepo_stats_ch,
        gene_loc_ch      // Pass gene coordinates for filtering
    )
    
    // Step 2: Run MAGMA analysis
    RUN_MAGMA(
        gwas_ch,
        gene_loc_ch,
        CREATE_MAGMA_GENESET.out.gene_set,
        magma_bin_ch,
        ref_files_ch,
        genome_build_ch
    )
    
    emit:
    gsa_results = RUN_MAGMA.out.gsa_file
    gene_results = RUN_MAGMA.out.genes_file
    genes_raw = RUN_MAGMA.out.genes_raw
    log_file = RUN_MAGMA.out.log_file
}