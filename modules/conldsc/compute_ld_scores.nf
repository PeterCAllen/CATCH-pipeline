// modules/conldsc/compute_ld_scores.nf
// CELLECT rules: compute_LD_scores, compute_LD_scores_all_genes

process COMPUTE_LD_SCORES {
    tag "${specificity_id}:chr${chr}"
    label 'medium_mem'

    container "${projectDir}/environments/ldsc-timshel.sif"

    input:
    tuple val(specificity_id), val(chr), path(annot), path(plink_files), path(print_snps)

    output:
    tuple val(specificity_id), val(chr),
          path("${specificity_id}.COMBINED_ANNOT.${chr}.l2.ldscore.gz"),
          path("${specificity_id}.COMBINED_ANNOT.${chr}.l2.M"),
          path("${specificity_id}.COMBINED_ANNOT.${chr}.l2.M_5_50"),
          path(annot), emit: ldscores
    path "*.log", emit: log

    script:
    def bfile = plink_files.find { it.name.endsWith('.bed') }.name.replaceAll(/\.bed$/, '')
    """
    ldsc.py \\
        --l2 \\
        --bfile "${bfile}" \\
        --ld-wind-cm ${params.conldsc_ld_wind_cm} \\
        --annot "${annot}" \\
        --thin-annot \\
        --out "${specificity_id}.COMBINED_ANNOT.${chr}" \\
        --print-snps "${print_snps}" \\
        2>&1 | tee "${specificity_id}.COMBINED_ANNOT.${chr}.l2.stdout.log"

    for ext in l2.ldscore.gz l2.M l2.M_5_50; do
        if [ ! -s "${specificity_id}.COMBINED_ANNOT.${chr}.\${ext}" ]; then
            echo "ERROR: ldsc.py --l2 did not produce ${specificity_id}.COMBINED_ANNOT.${chr}.\${ext}"
            exit 1
        fi
    done
    """
}

process COMPUTE_LD_SCORES_ALL_GENES {
    tag "${specificity_id}:chr${chr}"
    label 'medium_mem'
    publishDir { "${params.outdir}/conldsc/${specificity_id}/control" }, mode: 'copy'

    container "${projectDir}/environments/ldsc-timshel.sif"

    input:
    tuple val(specificity_id), val(chr), path(annot), path(plink_files), path(print_snps)

    output:
    tuple val(specificity_id),
          path("all_genes_in_${specificity_id}.${chr}.l2.ldscore.gz"),
          path("all_genes_in_${specificity_id}.${chr}.l2.M"),
          path("all_genes_in_${specificity_id}.${chr}.l2.M_5_50"),
          path(annot), emit: ldscores
    path "*.log", emit: log

    script:
    def bfile = plink_files.find { it.name.endsWith('.bed') }.name.replaceAll(/\.bed$/, '')
    """
    ldsc.py \\
        --l2 \\
        --bfile "${bfile}" \\
        --ld-wind-cm ${params.conldsc_ld_wind_cm} \\
        --annot "${annot}" \\
        --thin-annot \\
        --out "all_genes_in_${specificity_id}.${chr}" \\
        --print-snps "${print_snps}" \\
        2>&1 | tee "all_genes_in_${specificity_id}.${chr}.l2.stdout.log"

    for ext in l2.ldscore.gz l2.M l2.M_5_50; do
        if [ ! -s "all_genes_in_${specificity_id}.${chr}.\${ext}" ]; then
            echo "ERROR: ldsc.py --l2 did not produce all_genes_in_${specificity_id}.${chr}.\${ext}"
            exit 1
        fi
    done
    """
}
