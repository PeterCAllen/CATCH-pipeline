process RUN_SLDSC {
    tag "${cell_type}"
    label 'medium_mem'
    publishDir "${params.outdir}/ldsc/results", mode: 'copy'
    
    container "${projectDir}/environments/ldsc_v1.0.1.sif"

    input:
    tuple val(cell_type), 
          path(gwas_munged), 
          val(genome_build),
          path(baseline_dir), 
          path(weights_dir), 
          path(plink_dir),
          val(plink_prefix), 
          path(hapmap3), 
          path(annot_files)

    output:
    path "${cell_type}_sldsc.results", emit: results
    path "${cell_type}_sldsc.log", emit: log
    path "${cell_type}_sldsc.part_delete", emit: delete_vals, optional: true

    script:
    // Determine the correct prefix based on genome build
    def frq_prefix = (genome_build in ['hg38', 'GRCh38']) ? '1000G.EUR.hg38' : '1000G.EUR.hg19'
    
    """
    echo "========================================"
    echo "STEP 4: Stratified LD Score Regression"
    echo "========================================"
    echo "Cell type: ${cell_type}"
    echo "Genome build: ${genome_build}"
    echo "GWAS munged: ${gwas_munged}"
    echo ""

    # Stage baseline annotation files
    echo "Staging baseline annotations..."
    for chr in {1..22}; do
        for ext in annot.gz l2.ldscore.gz l2.M l2.M_5_50; do
            if [ -f "${baseline_dir}/baselineLD.\${chr}.\${ext}" ]; then
                ln -sf "${baseline_dir}/baselineLD.\${chr}.\${ext}" .
            fi
        done
    done

    # Stage weights files
    echo "Staging weights..."
    for chr in {1..22}; do
        for ext in l2.ldscore.gz; do
            if [ -f "${weights_dir}/weights.hm3_noMHC.\${chr}.\${ext}" ]; then
                ln -sf "${weights_dir}/weights.hm3_noMHC.\${chr}.\${ext}" .
            fi
        done
    done

    # Stage frequency files
    echo "Staging frequency files..."
    for chr in {1..22}; do
        if [ -f "${plink_dir}/${frq_prefix}.\${chr}.frq" ]; then
            ln -sf "${plink_dir}/${frq_prefix}.\${chr}.frq" .
        fi
    done

    # Annotation files are already staged in work directory
    echo ""
    echo "Annotation files present:"
    echo "  Annot files (.annot.gz): \$(ls ${cell_type}.*.annot.gz 2>/dev/null | wc -l)"
    echo "  LD score files (.l2.ldscore.gz): \$(ls ${cell_type}.*.l2.ldscore.gz 2>/dev/null | wc -l)"
    echo "  L2.M files (.l2.M): \$(ls ${cell_type}.*.l2.M 2>/dev/null | wc -l)"
    echo "  L2.M_5_50 files (.l2.M_5_50): \$(ls ${cell_type}.*.l2.M_5_50 2>/dev/null | wc -l)"
    echo ""

    # Verify all required files are present
    echo "Checking for required files..."
    MISSING=0
    for chr in {1..22}; do
        # Check baseline files
        if [ ! -f "baselineLD.\${chr}.l2.ldscore.gz" ]; then
            echo "WARNING: Missing baselineLD.\${chr}.l2.ldscore.gz"
            MISSING=1
        fi
        if [ ! -f "baselineLD.\${chr}.l2.M" ]; then
            echo "WARNING: Missing baselineLD.\${chr}.l2.M"
            MISSING=1
        fi
        if [ ! -f "baselineLD.\${chr}.l2.M_5_50" ]; then
            echo "WARNING: Missing baselineLD.\${chr}.l2.M_5_50"
            MISSING=1
        fi
        
        # Check weights
        if [ ! -f "weights.hm3_noMHC.\${chr}.l2.ldscore.gz" ]; then
            echo "WARNING: Missing weights.hm3_noMHC.\${chr}.l2.ldscore.gz"
            MISSING=1
        fi
        
        # Check frequency
        if [ ! -f "${frq_prefix}.\${chr}.frq" ]; then
            echo "WARNING: Missing ${frq_prefix}.\${chr}.frq"
            MISSING=1
        fi
        
        # Check annotation files
        if [ ! -f "${cell_type}.\${chr}.annot.gz" ]; then
            echo "WARNING: Missing ${cell_type}.\${chr}.annot.gz"
            MISSING=1
        fi
        if [ ! -f "${cell_type}.\${chr}.l2.ldscore.gz" ]; then
            echo "WARNING: Missing ${cell_type}.\${chr}.l2.ldscore.gz"
            MISSING=1
        fi
        if [ ! -f "${cell_type}.\${chr}.l2.M" ]; then
            echo "WARNING: Missing ${cell_type}.\${chr}.l2.M"
            MISSING=1
        fi
        if [ ! -f "${cell_type}.\${chr}.l2.M_5_50" ]; then
            echo "WARNING: Missing ${cell_type}.\${chr}.l2.M_5_50"
            MISSING=1
        fi
    done

    if [ \$MISSING -eq 1 ]; then
        echo "ERROR: Some required files are missing. Check warnings above."
        
    fi

    echo "All required files present. Running ldsc.py..."
    echo ""

    # Run LDSC
    ldsc.py \\
        --h2 ${gwas_munged} \\
        --ref-ld-chr baselineLD.,${cell_type}. \\
        --w-ld-chr weights.hm3_noMHC. \\
        --frqfile-chr ${frq_prefix}. \\
        --overlap-annot \\
        --print-coefficients \\
        --print-delete-vals \\
        --out ${cell_type}_sldsc \\
        2>&1 | tee ${cell_type}_sldsc.log

    echo ""
    echo "LDSC completed for ${cell_type}"
    """
}