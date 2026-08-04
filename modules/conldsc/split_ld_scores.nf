// modules/conldsc/split_ld_scores.nf
// CELLECT rule: split_LD_scores

process SPLIT_LD_SCORES {
    tag "${specificity_id}:chr${chr}"
    label 'medium_mem'

    container "${projectDir}/environments/cellect-py3.sif"

    input:
    tuple val(specificity_id), val(chr), path(ldscore), path(m), path(m_5_50), path(annot)

    output:
    tuple val(specificity_id), path("per_annotation/*"), emit: per_annotation

    script:
    """
    mkdir -p per_annotation

    python3 -u "${projectDir}/bin/cellect_split_ldscores.py" \\
        --ldscore "${ldscore}" \\
        --out_dir per_annotation
    """
}
