// modules/make_ldsc_continuous_annot.nf

process MAKE_LDSC_CONTINUOUS_ANNOT {
    tag "${cell_type}_chr${chr}"
    label 'process_medium'
    publishDir "${params.outdir}/ldsc/annotations", mode: 'copy', pattern: "*.{annot.gz,l2.ldscore.gz,l2.M,l2.M_5_50,log}"
    
    // Error handling: Continue pipeline even if one chromosome fails
    errorStrategy 'ignore'
    maxRetries 2
    
    container "file://${projectDir}/environments/py-r-cepo-scdrs.sif"

    input:
    tuple val(cell_type), val(chr), path(bed_file)
    val genome_build
    path plink_files  // All PLINK files for this chromosome (.bed, .bim, .fam)
    path hapmap3_file

    output:
    tuple val(cell_type), val(chr), path("*.l2.ldscore.gz"), emit: ldscore_file
    tuple val(cell_type), val(chr), path("*.annot.gz"), emit: annot_file
    tuple val(cell_type), val(chr), path("*.l2.M"), emit: l2_m, optional: true
    tuple val(cell_type), val(chr), path("*.l2.M_5_50"), emit: l2_m_5_50, optional: true
    path "*.log", emit: log

    script:
    def prefix = "${cell_type}.${chr}"
    
    // The plink_files are staged, so we use the basename
    def plink_prefix = (genome_build in ['hg38', 'GRCh38']) 
        ? "1000G.EUR.hg38"
        : "1000G.EUR.hg19"
    
    def bim_file = "${plink_prefix}.${chr}.bim"
    def bfile = "${plink_prefix}.${chr}"
    
    // LDSC parameters
    def ld_wind_cm = params.ldsc_ld_wind_cm ?: 1
    def thin_annot_flag = params.ldsc_thin_annot ? '--thin-annot' : ''
    
    """
    echo "========================================"
    echo "STEP 3: Create LDSC Annotation"
    echo "========================================"
    echo "Cell type: ${cell_type}"
    echo "Chromosome: ${chr}"
    echo "Genome build: ${genome_build}"
    echo "BED file: ${bed_file}"
    echo "BIM file: ${bim_file}"
    echo "PLINK bfile: ${bfile}"
    echo "HapMap3 SNPs: ${hapmap3_file}"
    echo ""
    
    # Verify input files exist
    if [ ! -f "${bed_file}" ]; then
        echo "❌ ERROR: BED file not found: ${bed_file}"
        exit 1
    fi
    
    if [ ! -f "${bim_file}" ]; then
        echo "❌ ERROR: BIM file not found: ${bim_file}"
        exit 1
    fi
    
    if [ ! -f "${hapmap3_file}" ]; then
        echo "❌ ERROR: HapMap3 file not found: ${hapmap3_file}"
        exit 1
    fi
    
    # Step 3a: Create annotation file using R script
    echo "Creating annotation file..."
    Rscript ${projectDir}/bin/ldsc-2-make_ldsc_continuous_annot.R \\
        "${bed_file}" \\
        "${bim_file}" \\
        "${prefix}.annot.gz" \\
        "full" \\
        "${cell_type}"
    
    # Verify annotation was created
    if [ ! -f "${prefix}.annot.gz" ]; then
        echo "❌ ERROR: Annotation file not created: ${prefix}.annot.gz"
        exit 1
    fi
    
    echo "✓ Annotation file created"
    echo ""
    
    # Step 3b: Compute LD scores
    echo "Computing LD scores with ldsc.py..."
    ldsc.py \\
        --l2 \\
        --bfile "${bfile}" \\
        --ld-wind-cm ${ld_wind_cm} \\
        --annot "${prefix}.annot.gz" \\
        --out "${prefix}" \\
        --print-snps "${hapmap3_file}" \\
        ${thin_annot_flag} \\
        2>&1 | tee ${prefix}.log
    
    # Verify LD scores were created
    if [ ! -f "${prefix}.l2.ldscore.gz" ]; then
        echo "❌ ERROR: LD score file not created: ${prefix}.l2.ldscore.gz"
        echo "Check ${prefix}.log for details"
        exit 1
    fi
    
    echo ""
    echo "✓ LD scores computed successfully"
    echo "Files created:"
    ls -lh ${prefix}.*
    """
}
