// modules/create_magma_geneset.nf

process CREATE_MAGMA_GENESET {
    tag "$cepo_stats"
    container "${projectDir}/environments/py-r-cepo-scdrs.sif"
    publishDir "${params.outdir}/magma", mode: 'copy'

    input:
    path cepo_stats
    path magma_gene_loc  // MAGMA gene location file (Gene, chr, start, end, strand, gene_name)

    output:
    path "magma.sets.txt", emit: gene_set

    script:
    """
    Rscript "${projectDir}/bin/cepo_to_magma_geneset.R" \\
        --cepo_stats "${cepo_stats}" \\
        --magma_gene_loc "${magma_gene_loc}" \\
        --top_pct ${params.magma_top_pct ?: 0.10} \\
        --output_file magma.sets.txt
    """
}
