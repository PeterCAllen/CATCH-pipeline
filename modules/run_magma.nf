process RUN_MAGMA {
    tag "$gwas"
    container "${projectDir}/environments/py-r-cepo-scdrs.sif"
    publishDir "${params.outdir}/magma", mode: 'copy', pattern: "*.{gsa.out,genes.out,genes.raw,log}"

    input:
    path gwas
    path gene_loc
    path gene_set
    path magma_bin
    path "*"  // Stage all reference files (chr 1-22 .bed, .bim, .fam) directly into work directory
    val genome_build

    output:
    path "*.gsa.out", emit: gsa_file
    path "*.genes.out", emit: genes_file
    path "*.genes.raw", emit: genes_raw
    path "*.log", emit: log_file

    script:
    def gwas_prefix = gwas.getBaseName().replaceAll(/\.(txt\.gz|txt|gz)$/, '')
    
    // Reference prefix for per-chromosome files (MAGMA will find them automatically)
    // Files are staged directly in work directory, we'll use them from there
    def ref_prefix = (genome_build in ['hg38', 'GRCh38']) ? 
        '1000G.EUR.hg38' : 
        '1000G.EUR.QC'
    
    """
    echo "================================================================"
    echo "MAGMA Gene-Set Analysis Pipeline"
    echo "================================================================"
    echo "GWAS file: ${gwas}"
    echo "Genome build: ${genome_build}"
    echo "Sample size: ${params.gwas_sample_size}"
    echo "================================================================"
    echo ""

    # Input validation
    echo "Checking required inputs..."
    echo "----------------------------"
    
    if [ ! -f "${gwas}" ]; then
        echo "❌ ERROR: GWAS summary statistics not found: ${gwas}"
        exit 1
    fi
    echo "✓ Found GWAS file: ${gwas}"
    GWAS_LINES=\$(zcat -f "${gwas}" | wc -l)
    echo "  - Lines: \${GWAS_LINES}"
    
    if [ ! -f "${gene_loc}" ]; then
        echo "❌ ERROR: MAGMA gene coordinates not found: ${gene_loc}"
        exit 1
    fi
    echo "✓ Found MAGMA coordinates: ${gene_loc}"
    
    if [ ! -f "${gene_set}" ]; then
        echo "❌ ERROR: MAGMA gene sets not found: ${gene_set}"
        exit 1
    fi
    echo "✓ Found MAGMA gene sets: ${gene_set}"
    NSETS=\$(wc -l < "${gene_set}")
    echo "  - Gene sets: \${NSETS}"
    
    # Check for per-chromosome reference files
    if [ ! -f "${ref_prefix}.1.bed" ]; then
        echo "❌ ERROR: Reference data not found: ${ref_prefix}.*.bed"
        echo "Looking for files in work directory:"
        ls -la *.bed | head -n 10 || true
        exit 1
    fi
    echo "✓ Found per-chromosome reference data: ${ref_prefix}.*.{bed,bim,fam}"
    
    if [ ! -f "./${magma_bin}" ]; then
        echo "❌ ERROR: MAGMA binary not found: ${magma_bin}"
        exit 1
    fi
    echo "✓ Found MAGMA binary: ${magma_bin}"
    chmod +x "./${magma_bin}"

    echo ""
    echo "========================================"
    echo "Step 1: Format GWAS for MAGMA"
    echo "========================================"
    echo ""
    
    # The combined BIM file is staged into work directory - find it
    # It will have the full name like "1000G.EUR.hg38.bim" or "1000G.EUR.QC.bim"
    COMBINED_BIM="\$(ls ${ref_prefix}.bim 2>/dev/null || echo "")"
    
    if [ -z "\$COMBINED_BIM" ] || [ ! -f "\$COMBINED_BIM" ]; then
        echo "❌ ERROR: Combined BIM file not found: ${ref_prefix}.bim"
        echo "Files in work directory:"
        ls -la *.bim | head -n 10 || true
        exit 1
    fi
    
    BIM_SNPS=\$(wc -l < "\$COMBINED_BIM")
    echo "✓ Found combined BIM file: \$COMBINED_BIM with \${BIM_SNPS} SNPs"

    # Check if GWAS file has header
    FIRST_LINE=\$(zcat -f "${gwas}" | head -n 1)
    if echo "\$FIRST_LINE" | grep -Eq '^(SNP|rsid|MarkerName|SNPID|CHR)'; then
        echo "Detected header in GWAS file"
    else
        echo "No header detected in GWAS file"
    fi

    # Format GWAS file using the formatting script
    echo "Formatting GWAS data for MAGMA..."
    bash ${projectDir}/bin/format_gwas_for_magma.sh \
        "${gwas}" \
        "\$COMBINED_BIM" \
        "${gwas_prefix}.formatted"

    # Validate formatted outputs
    if [ ! -f "${gwas_prefix}.formatted.pval" ]; then
        echo "❌ ERROR: Failed to format GWAS file (missing .pval)"
        exit 1
    fi
    
    if [ ! -f "${gwas_prefix}.formatted.snp.loc" ]; then
        echo "❌ ERROR: Failed to format GWAS file (missing .snp.loc)"
        exit 1
    fi
    
    echo "✓ Created formatted GWAS files"
    echo "Preview of .pval file:"
    head -n 5 "${gwas_prefix}.formatted.pval"
    echo ""
    echo "Preview of .snp.loc file:"
    head -n 5 "${gwas_prefix}.formatted.snp.loc"

    echo ""
    echo "========================================"
    echo "Step 2: Annotate SNPs to Genes"
    echo "========================================"
    echo ""
    
    echo "Running MAGMA annotation..."
    echo "This may take several minutes..."
    
    ./${magma_bin} \
        --annotate window=${params.magma_window_kb ?: 10} \
        --snp-loc "${gwas_prefix}.formatted.snp.loc" \
        --gene-loc "${gene_loc}" \
        --out "${gwas_prefix}_annot"
    
    if [ ! -f "${gwas_prefix}_annot.genes.annot" ]; then
        echo "❌ ERROR: MAGMA annotation failed"
        exit 1
    fi
    
    echo "✓ SNP annotation complete"
    ANNOT_GENES=\$(grep -v "^#" "${gwas_prefix}_annot.genes.annot" | wc -l)
    echo "  - Genes with SNPs: \${ANNOT_GENES}"

    echo ""
    echo "========================================"
    echo "Step 3: Gene Analysis"
    echo "========================================"
    echo ""
    
    echo "Running MAGMA gene analysis..."
    echo "This may take several minutes..."
    
    echo "Preparing reference files for MAGMA..."
    # MAGMA expects: --bfile prefix merge
    # where prefix is the common prefix and MAGMA will look for prefix.1.bed, prefix.2.bed, etc.
    # Our files are named: 1000G.EUR.hg38.1.bed, 1000G.EUR.hg38.2.bed, etc.
    # So we use the prefix without the chromosome number
    
    echo "Using reference prefix: ${ref_prefix}"
    echo "MAGMA will look for files: ${ref_prefix}.{1..22}.{bed,bim,fam}"
    
    # Verify first chromosome file exists
    if [ ! -f "${ref_prefix}.1.bed" ]; then
        echo "❌ ERROR: Cannot find ${ref_prefix}.1.bed"
        echo "Files in work directory:"
        ls -la *.bed | head -n 20 || true
        exit 1
    fi
    echo "✓ Verified chromosome 1 reference files exist"
    
    ./${magma_bin} \
        --bfile ${ref_prefix} \
        --gene-annot "${gwas_prefix}_annot.genes.annot" \
        --pval "${gwas_prefix}.formatted.pval" ncol=N \
        --out "${gwas_prefix}.magma_genes"
    
    if [ ! -f "${gwas_prefix}.magma_genes.genes.out" ]; then
        echo "❌ ERROR: MAGMA gene analysis failed"
        exit 1
    fi
    
    echo "✓ Gene analysis complete"
    echo "Preview of top genes:"
    head -n 10 "${gwas_prefix}.magma_genes.genes.out"
    
    SIG_GENES=\$(awk '\$9 < 0.05' "${gwas_prefix}.magma_genes.genes.out" | wc -l)
    echo ""
    echo "Significant genes (P < 0.05): \${SIG_GENES}"

    echo ""
    echo "========================================"
    echo "Step 4: Gene-Set Analysis"
    echo "========================================"
    echo ""
    
    echo "Running MAGMA gene-set analysis..."
    
    ./${magma_bin} \
        --gene-results "${gwas_prefix}.magma_genes.genes.raw" \
        --set-annot "${gene_set}" \
        --out "${gwas_prefix}.magma_geneset"
    
    if [ ! -f "${gwas_prefix}.magma_geneset.gsa.out" ]; then
        echo "❌ ERROR: MAGMA gene-set analysis failed"
        exit 1
    fi
    
    echo "✓ Gene-set analysis complete"
    echo ""
    echo "Results:"
    head -n 10 "${gwas_prefix}.magma_geneset.gsa.out"

    echo ""
    echo "========================================"
    echo "RESULTS SUMMARY"
    echo "========================================"
    echo ""
    
    SIG_SETS=\$(awk 'NR>1 && \$7 < 0.05' "${gwas_prefix}.magma_geneset.gsa.out" | wc -l)
    echo "Significant gene sets (P < 0.05): \${SIG_SETS} / \${NSETS}"
    
    echo ""
    echo "Output files created:"
    ls -lh "${gwas_prefix}"*.{out,annot,raw,pval} 2>/dev/null || true

    echo ""
    echo "========================================"
    echo "VALIDATION"
    echo "========================================"
    echo ""
    
    EXPECTED_FILES=(
        "${gwas_prefix}.formatted.pval"
        "${gwas_prefix}.formatted.snp.loc"
        "${gwas_prefix}_annot.genes.annot"
        "${gwas_prefix}.magma_genes.genes.out"
        "${gwas_prefix}.magma_genes.genes.raw"
        "${gwas_prefix}.magma_geneset.gsa.out"
    )
    
    ALL_EXIST=true
    for f in "\${EXPECTED_FILES[@]}"; do
        if [ ! -f "\$f" ]; then
            echo "❌ Missing: \$(basename \$f)"
            ALL_EXIST=false
        else
            echo "✓ Found: \$(basename \$f)"
        fi
    done
    
    echo ""
    if [ "\$ALL_EXIST" = true ]; then
        echo "✓ All expected MAGMA output files created successfully!"
    else
        echo "⚠️  Some output files are missing"
        exit 1
    fi

    echo ""
    echo "========================================"
    echo "MAGMA Complete!"
    echo "========================================"
    echo ""
    echo "Summary:"
    echo "  - Genome build: ${genome_build}"
    echo "  - GWAS SNPs processed: \${GWAS_LINES}"
    echo "  - Genes annotated: \${ANNOT_GENES}"
    echo "  - Significant genes: \${SIG_GENES}"
    echo "  - Gene sets tested: \${NSETS}"
    echo "  - Significant gene sets: \${SIG_SETS}"
    echo ""
    echo "Next: Review results in ${gwas_prefix}.magma_geneset.gsa.out for cell type enrichments in your GWAS trait"
    echo "================================================================"
    """
}