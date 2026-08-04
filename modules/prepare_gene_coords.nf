// modules/prepare_gene_coords.nf

include { asBool } from '../lib/util'

process PREPARE_GENE_COORDS {
    tag "${genome_build}"
    publishDir "${params.outdir}/reference", mode: 'copy'
    container "${projectDir}/environments/py-r-cepo-scdrs.sif"

    input:
    path gene_matrix
    val genome_build

    output:
    path "magma_coords.${genome_build}.loc",   emit: magma_loc
    path "ldsc_coords.${genome_build}.loc",    emit: ldsc_loc
    path "cepo_coords.${genome_build}.tsv",    emit: cepo_coords
    path "cellect_coords.${genome_build}.loc", emit: cellect_loc

    script:
    def pc_flag = asBool(params.protein_coding_only) ? '' : '--no_protein_coding'
    """
    python3 -u "${projectDir}/bin/generate_gene_coords.py" \\
        --gene_matrix "${gene_matrix}" \\
        --genome_build "${genome_build}" \\
        --output_magma "magma_coords.${genome_build}.loc" \\
        --output_ldsc "ldsc_coords.${genome_build}.loc" \\
        --output_cepo "cepo_coords.${genome_build}.tsv" \\
        --output_cellect "cellect_coords.${genome_build}.loc" \\
        ${pc_flag}
    """
}
