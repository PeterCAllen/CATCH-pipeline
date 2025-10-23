// modules/prepare_mbat_genes.nf
// Convert gene coordinates to mBAT format

process PREPARE_MBAT_GENES {
    tag "mbat_gene_list"
    label 'low_mem'
    publishDir "${params.outdir}/scdrs/gwas", mode: 'copy'
    
    input:
    path gene_coords

    output:
    path "genes.mbat.loc", emit: mbat_genes
    path "*.log", emit: log

    script:
    """
    echo "Converting gene coordinates to mBAT format..." | tee prepare_genes.log
    echo "Input: ${gene_coords}" | tee -a prepare_genes.log
    echo "" | tee -a prepare_genes.log
    
    # mBAT needs: Gene, Chr, Start, End (no header, tab-separated)
    tail -n +2 "${gene_coords}" | awk -F'\\t' '{print \$2 "\\t" \$3 "\\t" \$4 "\\t" \$1}' > genes.mbat.loc
    
    N_GENES=\$(wc -l < genes.mbat.loc)
    echo "✓ Created mBAT gene list with \${N_GENES} genes" | tee -a prepare_genes.log
    echo "" | tee -a prepare_genes.log
    echo "First 5 genes:" | tee -a prepare_genes.log
    head -n 5 genes.mbat.loc | tee -a prepare_genes.log
    
    mv prepare_genes.log genes.mbat.log
    """
}
