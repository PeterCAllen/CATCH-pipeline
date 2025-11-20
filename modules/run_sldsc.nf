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
          path(frq_dir),           // FIXED: Changed from plink_dir to frq_dir
          val(plink_prefix),       // This is just the prefix name (e.g., "1000G.EUR.QC")
          path(hapmap3), 
          path(annot_files)

    output:
    path "${cell_type}_sldsc.results", emit: results
    path "${cell_type}_sldsc.log", emit: log
    path "${cell_type}_sldsc.part_delete", emit: delete_vals, optional: true

    script:
    // Extract the baseline prefix name from params (e.g., "baseline." from "data/.../baseline.")
    def baseline_prefix = params.ref_hg19_baseline.substring(params.ref_hg19_baseline.lastIndexOf('/') + 1)
    // Remove trailing dot if present for use in filenames
    def baseline_name = baseline_prefix.endsWith('.') ? baseline_prefix.substring(0, baseline_prefix.length() - 1) : baseline_prefix
    
    // Use the plink_prefix value that was passed in (already just the filename)
    def frq_prefix = plink_prefix
    
    """
    echo "========================================"
    echo "STEP 4: Stratified LD Score Regression"
    echo "========================================"
    echo "Cell type: ${cell_type}"
    echo "Genome build: ${genome_build}"
    echo "GWAS munged: ${gwas_munged}"
    echo "Baseline prefix: ${baseline_name}"
    echo "Frequency prefix: ${frq_prefix}"
    echo ""

    # Stage baseline annotation files
    echo "Staging baseline annotations..."
    for chr in {1..22}; do
        for ext in annot.gz l2.ldscore.gz l2.M l2.M_5_50; do
            if [ -f "${baseline_dir}/${baseline_name}.\${chr}.\${ext}" ]; then
                ln -sf "${baseline_dir}/${baseline_name}.\${chr}.\${ext}" .
            else
                echo "WARNING: Missing ${baseline_dir}/${baseline_name}.\${chr}.\${ext}"
            fi
        done
    done

    # Stage weights files
    echo "Staging weights..."
    for chr in {1..22}; do
        for ext in l2.ldscore.gz; do
            if [ -f "${weights_dir}/weights.hm3_noMHC.\${chr}.\${ext}" ]; then
                ln -sf "${weights_dir}/weights.hm3_noMHC.\${chr}.\${ext}" .
            else
                echo "WARNING: Missing ${weights_dir}/weights.hm3_noMHC.\${chr}.\${ext}"
            fi
        done
    done

    # Stage frequency files from frq_dir
    echo "Staging frequency files from ${frq_dir}..."
    for chr in {1..22}; do
        if [ -f "${frq_dir}/${frq_prefix}.\${chr}.frq" ]; then
            ln -sf "${frq_dir}/${frq_prefix}.\${chr}.frq" .
        else
            echo "WARNING: Missing ${frq_dir}/${frq_prefix}.\${chr}.frq"
        fi
    done

    # Verify staged files
    echo ""
    echo "Staged files in working directory:"
    echo "  Baseline files: \$(ls ${baseline_name}.*.l2.ldscore.gz 2>/dev/null | wc -l)"
    echo "  Weight files: \$(ls weights.hm3_noMHC.*.l2.ldscore.gz 2>/dev/null | wc -l)"
    echo "  Frequency files: \$(ls ${frq_prefix}.*.frq 2>/dev/null | wc -l)"
    echo ""

    # Annotation files are already staged in work directory
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
        if [ ! -f "${baseline_name}.\${chr}.l2.ldscore.gz" ]; then
            echo "ERROR: Missing ${baseline_name}.\${chr}.l2.ldscore.gz"
            MISSING=1
        fi
        if [ ! -f "${baseline_name}.\${chr}.l2.M" ]; then
            echo "ERROR: Missing ${baseline_name}.\${chr}.l2.M"
            MISSING=1
        fi
        if [ ! -f "${baseline_name}.\${chr}.l2.M_5_50" ]; then
            echo "ERROR: Missing ${baseline_name}.\${chr}.l2.M_5_50"
            MISSING=1
        fi
        
        # Check weights
        if [ ! -f "weights.hm3_noMHC.\${chr}.l2.ldscore.gz" ]; then
            echo "ERROR: Missing weights.hm3_noMHC.\${chr}.l2.ldscore.gz"
            MISSING=1
        fi
        
        # Check frequency
        if [ ! -f "${frq_prefix}.\${chr}.frq" ]; then
            echo "ERROR: Missing ${frq_prefix}.\${chr}.frq"
            MISSING=1
        fi
        
        # Check annotation files
        if [ ! -f "${cell_type}.\${chr}.annot.gz" ]; then
            echo "ERROR: Missing ${cell_type}.\${chr}.annot.gz"
            MISSING=1
        fi
        if [ ! -f "${cell_type}.\${chr}.l2.ldscore.gz" ]; then
            echo "ERROR: Missing ${cell_type}.\${chr}.l2.ldscore.gz"
            MISSING=1
        fi
        if [ ! -f "${cell_type}.\${chr}.l2.M" ]; then
            echo "ERROR: Missing ${cell_type}.\${chr}.l2.M"
            MISSING=1
        fi
        if [ ! -f "${cell_type}.\${chr}.l2.M_5_50" ]; then
            echo "ERROR: Missing ${cell_type}.\${chr}.l2.M_5_50"
            MISSING=1
        fi
    done

    if [ \$MISSING -eq 1 ]; then
        echo "ERROR: Some required files are missing. Check errors above."
        exit 1
    fi

    echo "All required files present. Running ldsc.py..."
    echo ""

    # Run LDSC
    ldsc.py \\
        --h2 ${gwas_munged} \\
        --ref-ld-chr ${baseline_name}.,${cell_type}. \\
        --w-ld-chr weights.hm3_noMHC. \\
        --frqfile-chr ${frq_prefix}. \\
        --overlap-annot \\
        --print-coefficients \\
        --print-delete-vals \\
        --out ${cell_type}_sldsc \\
        2>&1 | tee ${cell_type}_sldsc.log

    if [ \$? -ne 0 ]; then
        echo "ERROR: ldsc.py failed"
        exit 1
    fi

    echo ""
    echo "LDSC completed for ${cell_type}"
    echo "Output files:"
    ls -lh ${cell_type}_sldsc.*
    """
}
