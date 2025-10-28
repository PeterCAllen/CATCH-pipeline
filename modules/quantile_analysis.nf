process COMPUTE_QUANTILE_M {
    tag "${cell_type}"
    label 'medium_mem'
    publishDir "${params.outdir}/ldsc/quantile_results", mode: 'copy', pattern: "*.q5.M"
    
    container "${projectDir}/environments/perl_latest.sif"

    input:
    tuple val(cell_type),
          path(annot_files),
          val(genome_build),
          path(baseline_dir),
          path(plink_dir),
          val(plink_prefix)

    output:
    tuple val(cell_type), path("${cell_type}.q5.M"), emit: quantile_m
    path "${cell_type}_quantile_m.log", emit: log

    script:
    def frq_prefix = (genome_build in ['hg38', 'GRCh38']) ? '1000G.EUR.hg38' : '1000G.EUR.hg19'
    
    """
    echo "========================================"
    echo "STEP 5a: Compute Quantile M Values"
    echo "========================================"
    echo "Cell type: ${cell_type}"
    echo "Genome build: ${genome_build}"
    echo ""

    # Stage baseline annotation files
    echo "Staging baseline annotations..."
    for chr in {1..22}; do
        if [ -f "${baseline_dir}/baselineLD.\${chr}.annot.gz" ]; then
            ln -sf "${baseline_dir}/baselineLD.\${chr}.annot.gz" .
        fi
    done

    # Stage frequency files
    echo "Staging frequency files..."
    for chr in {1..22}; do
        if [ -f "${plink_dir}/${frq_prefix}.\${chr}.frq" ]; then
            ln -sf "${plink_dir}/${frq_prefix}.\${chr}.frq" .
        fi
    done

    # Verify annotation files are present
    echo ""
    echo "Cell-type annotation files present:"
    ls -lh ${cell_type}.*.annot.gz 2>/dev/null || echo "WARNING: No annotation files found"
    echo ""

    # Check if Perl script exists
    if [ ! -f "${projectDir}/bin/ldsc-quantile-M.pl" ]; then
        echo "ERROR: Perl script not found at ${projectDir}/bin/ldsc-quantile-M.pl"
        
    fi

    echo "Running Perl script for quantile M calculation..."
    perl ${projectDir}/bin/ldsc-quantile-M.pl \\
        --ref-annot-chr baselineLD.,${cell_type}. \\
        --frqfile-chr ${frq_prefix}. \\
        --annot-header "ANNOT" \\
        --nb-quantile 5 \\
        --maf 0.05 \\
        --out ${cell_type}.q5.M \\
        2>&1 | tee ${cell_type}_quantile_m.log

    if [ \$? -ne 0 ]; then
        echo "ERROR: Perl script failed"
        
    fi

    echo ""
    echo "Quantile M calculation completed for ${cell_type}"
    echo "Output file:"
    ls -lh ${cell_type}.q5.M
    """
}

process COMPUTE_QUANTILE_H2G {
    tag "${cell_type}"
    label 'medium_mem'
    publishDir "${params.outdir}/ldsc/quantile_results", mode: 'copy'
    
    container "${projectDir}/environments/py-r-cepo-scdrs.sif"

    input:
    tuple val(cell_type),
          path(quantile_m),
          path(sldsc_results),
          path(sldsc_parts_delete)

    output:
    path "${cell_type}.q5.txt", emit: quantile_results
    path "${cell_type}_quantile_h2g.log", emit: log

    script:
    // Remove .results extension from sldsc_results filename
    def sldsc_base = sldsc_results.baseName
    
    """
    echo "========================================"
    echo "STEP 5b: Compute Quantile h2g"
    echo "========================================"
    echo "Cell type: ${cell_type}"
    echo "Quantile M file: ${quantile_m}"
    echo "S-LDSC results file: ${sldsc_results}"
    echo "S-LDSC results base: ${sldsc_base}"
    echo ""

    # Check if R script exists
    if [ ! -f "${projectDir}/bin/ldsc-quantile_h2g.R" ]; then
        echo "ERROR: R script not found at ${projectDir}/bin/ldsc-quantile_h2g.R"
        
    fi

    echo "Running R script for quantile h2g analysis..."
    echo "Note: Passing base name '${sldsc_base}' to R script (it will add .results extension)"
    Rscript ${projectDir}/bin/ldsc-quantile_h2g.R \\
        ${quantile_m} \\
        ${sldsc_base} \\
        ${cell_type}.q5.txt \\
        2>&1 | tee ${cell_type}_quantile_h2g.log

    if [ \$? -ne 0 ]; then
        echo "ERROR: R script failed"
        
    fi

    echo ""
    echo "Quantile h2g analysis completed for ${cell_type}"
    echo "Output file:"
    ls -lh ${cell_type}.q5.txt
    """
}

process COMPARE_ANNOTATIONS {
    label 'medium_mem'
    publishDir "${params.outdir}/ldsc/quantile_comparison", mode: 'copy'
    
    container "${projectDir}/environments/py-r-cepo-scdrs.sif"

    input:
    path(quantile_results)
    val(gwas_name)

    output:
    path "${gwas_name}_annotation_comparison.csv", emit: comparison_table
    path "${gwas_name}_enrichment_comparison.pdf", emit: enrichment_plot
    path "${gwas_name}_prop_h2g_comparison.pdf", emit: prop_h2g_plot
    path "${gwas_name}_enrichment_vs_prop.pdf", emit: scatter_plot
    path "${gwas_name}_quintile_heatmap.pdf", emit: heatmap, optional: true
    path "compare_annotations.log", emit: log

    script:
    """
    echo "========================================"
    echo "STEP 6: Compare Annotations"
    echo "========================================"
    echo "GWAS name: ${gwas_name}"
    echo ""
    
    # List all input files
    echo "Input quantile result files:"
    ls -lh *.q5.txt
    echo ""

    # Check if R script exists
    if [ ! -f "${projectDir}/bin/ldsc-compare_annotations.R" ]; then
        echo "ERROR: R script not found at ${projectDir}/bin/ldsc-compare_annotations.R"
        
    fi

    # Run comparison script
    # Results are in current directory, output to current directory
    Rscript ${projectDir}/bin/ldsc-compare_annotations.R \\
        . \\
        ${gwas_name} \\
        2>&1 | tee compare_annotations.log
    """
}