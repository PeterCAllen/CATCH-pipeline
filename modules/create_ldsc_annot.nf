// modules/create_ldsc_annot.nf

process CREATE_LDSC_ANNOT {
    tag "${cell_type}_chr${chr}"
    label 'medium_mem'
    publishDir "${params.outdir}/ldsc/annotations", mode: 'copy'
    
    // Error handling: Continue pipeline even if one chromosome fails
    // cache false
    errorStrategy 'ignore'
    maxRetries 2
    
    container "${projectDir}/environments/py-r-cepo-scdrs.sif"

    input:
    tuple val(cell_type), val(chr), path(bed_file)
    val genome_build
    path plink_bim  // Just the .bim file for this chromosome

    output:
    tuple val(cell_type), val(chr), path("*.annot.gz"), emit: annot_file

    script:
    def prefix = "${cell_type}.${chr}"
    // Extract PLINK prefix from params
    def plink_prefix = new File(params.ref_hg19_plink_prefix).name
    def bim_file = "${plink_prefix}.${chr}.bim"
    
    """
    echo "========================================"
    echo "STEP 3a: Create LDSC Annotation"
    echo "========================================"
    echo "Cell type: ${cell_type}"
    echo "Chromosome: ${chr}"
    echo "Genome build: ${genome_build}"
    echo "BED file: ${bed_file}"
    echo "BIM file: ${bim_file}"
    echo ""
    
    # Verify input files exist
    if [ ! -f "${bed_file}" ]; then
        echo "❌ ERROR: BED file not found: ${bed_file}"
        
    fi
    
    if [ ! -f "${bim_file}" ]; then
        echo "❌ ERROR: BIM file not found: ${bim_file}"
        
    fi
    
    # Create annotation file using R script
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
        
    fi
    
    echo "✓ Annotation file created"
    ls -lh ${prefix}.annot.gz
    """
}
