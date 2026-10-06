// modules/conldsc/format_genes.nf
// CELLECT rule: format_genes

process FORMAT_GENES {
    tag "${params.conldsc_window_kb}kb"
    label 'low_mem'
    publishDir "${params.outdir}/conldsc/shared/bed", mode: 'copy'

    container "${projectDir}/environments/cellect-py3.sif"

    input:
    path gene_coords
    path chr_sizes

    output:
    path "genes_plus_${params.conldsc_window_kb}kb.*.bed", emit: gene_beds
    path "*.log",                                          emit: log

    script:
    """
    python3 -u "${projectDir}/bin/cellect_format_genes.py" \\
        --gene_coords "${gene_coords}" \\
        --chr_sizes "${chr_sizes}" \\
        --out_dir . \\
        --windowsize_kb ${params.conldsc_window_kb} \\
        2>&1 | tee format_genes.log

    N_BED=\$(ls genes_plus_${params.conldsc_window_kb}kb.*.bed 2>/dev/null | wc -l)
    if [ "\$N_BED" -ne 22 ]; then
        echo "ERROR: expected 22 per-chromosome BED files, produced \$N_BED"
        exit 1
    fi
    """
}
