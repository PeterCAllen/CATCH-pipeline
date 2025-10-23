// modules/create_continuous_beds.nf

process CREATE_CONTINUOUS_BEDS {
    tag "${cepo_stats.simpleName}"
    label 'medium_mem'
    publishDir "${params.outdir}/ldsc/bedfiles", mode: 'copy'
    
    container "${projectDir}/environments/py-r-cepo-scdrs.sif"

    input:
    path cepo_stats
    path gene_coords
    val window_kb

    output:
    path "*.bed", emit: bed_dir

    script:
    """
    echo "========================================"
    echo "STEP 2: Convert CEPO Results to BED Files"
    echo "========================================"
    echo "CEPO stats: ${cepo_stats}"
    echo "Gene coordinates: ${gene_coords}"
    echo "Window size: ${window_kb} kb"
    echo ""
    
    # Run the BED conversion script
    Rscript ${projectDir}/bin/ldsc-1-csv_matrix_to_bed.R \\
        "${cepo_stats}" \\
        "${gene_coords}" \\
        ${window_kb} \\
        .
    
    echo ""
    echo "✓ BED files created:"
    ls -lh *.bed | head -n 5
    echo ""
    echo "Total BED files: \$(ls *.bed | wc -l)"
    """
}
