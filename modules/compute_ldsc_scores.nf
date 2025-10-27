// modules/compute_ldsc_scores.nf

process COMPUTE_LDSC_SCORES {
    tag "${cell_type}_chr${chr}"
    label 'medium_mem'
    publishDir "${params.outdir}/ldsc/annotations", mode: 'copy'
    
    container "${projectDir}/environments/ldsc_v1.0.1.sif"

    input:
    tuple val(cell_type), val(chr), path(annot_file)
    val genome_build
    path plink_files  // All PLINK files for this chromosome (.bed, .bim, .fam)
    path hapmap3_file

    output:
    tuple val(cell_type), val(chr), path("*.l2.ldscore.gz"), emit: ldscore_file
    tuple val(cell_type), val(chr), path("*.l2.M"), emit: l2_m, optional: true
    tuple val(cell_type), val(chr), path("*.l2.M_5_50"), emit: l2_m_5_50, optional: true
    path "*.log", emit: log

    script:
    def prefix = "${cell_type}.${chr}"
    def plink_prefix = (genome_build in ['hg38', 'GRCh38']) 
        ? "1000G.EUR.hg38"
        : "1000G.EUR.hg19"
    def bfile = "${plink_prefix}.${chr}"
    
    // LDSC parameters
    def ld_wind_cm = params.ldsc_ld_wind_cm ?: 1
    
    """
    echo "========================================"
    echo "STEP 3b: Compute LD Scores"
    echo "========================================"
    echo "Cell type: ${cell_type}"
    echo "Chromosome: ${chr}"
    echo "Genome build: ${genome_build}"
    echo "Annotation file: ${annot_file}"
    echo "PLINK bfile: ${bfile}"
    echo "HapMap3 SNPs: ${hapmap3_file}"
    echo ""
    
    # Verify input files exist
    if [ ! -f "${annot_file}" ]; then
        echo "❌ ERROR: Annotation file not found: ${annot_file}"
        exit 1
    fi
    
    if [ ! -f "${bfile}.bim" ]; then
        echo "❌ ERROR: PLINK BIM file not found: ${bfile}.bim"
        exit 1
    fi
    
    if [ ! -f "${hapmap3_file}" ]; then
        echo "❌ ERROR: HapMap3 file not found: ${hapmap3_file}"
        exit 1
    fi
    
    # Compute LD scores
    echo "Computing LD scores with ldsc.py..."
    ldsc.py \\
        --l2 \\
        --bfile "${bfile}" \\
        --ld-wind-cm ${ld_wind_cm} \\
        --annot "${annot_file}" \\
        --out "${prefix}" \\
        --print-snps "${hapmap3_file}" \\
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
