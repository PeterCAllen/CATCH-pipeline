// modules/conldsc/make_annot.nf
// CELLECT rules: make_annot, make_annot_all_genes

include { asBool } from '../../lib/util'

process MAKE_ANNOT {
    tag "${specificity_id}:chr${chr}"
    label 'medium_mem'
    publishDir { "${params.outdir}/conldsc/${specificity_id}/annots" }, mode: 'copy', pattern: "*.keep.annot.gz"

    container "${projectDir}/environments/cellect-py3.sif"

    input:
    tuple val(specificity_id), val(chr), path(spec_matrix), path(overlap_segments), path(bim)

    output:
    tuple val(specificity_id), val(chr), path("${specificity_id}.COMBINED_ANNOT.${chr}.annot.gz"), emit: annot
    path "${specificity_id}.COMBINED_ANNOT.${chr}.keep.annot.gz", emit: annot_keep, optional: true

    script:
    def keep = asBool(params.conldsc_keep_annots) ? "--keep_annots --out_annot_keep ${specificity_id}.COMBINED_ANNOT.${chr}.keep.annot.gz" : ''
    """
    export TMPDIR=\${TMPDIR:-\$PWD/pybedtools_tmp}
    mkdir -p "\$TMPDIR"

    python3 -u "${projectDir}/bin/cellect_make_annot.py" \\
        --spec_matrix "${spec_matrix}" \\
        --overlap_segments "${overlap_segments}" \\
        --bimfile "${bim}" \\
        --chromosome ${chr} \\
        --out_dir . \\
        --out_prefix "${specificity_id}" \\
        ${keep}
    """
}

process MAKE_ANNOT_ALL_GENES {
    tag "${specificity_id}:chr${chr}"
    label 'medium_mem'

    container "${projectDir}/environments/cellect-py3.sif"

    input:
    tuple val(specificity_id), val(chr), path(all_genes), path(overlap_segments), path(bim)

    output:
    tuple val(specificity_id), val(chr), path("all_genes_in_${specificity_id}.${chr}.annot.gz"), emit: annot

    script:
    """
    export TMPDIR=\${TMPDIR:-\$PWD/pybedtools_tmp}
    mkdir -p "\$TMPDIR"

    python3 -u "${projectDir}/bin/cellect_make_annot.py" \\
        --spec_matrix "${all_genes}" \\
        --overlap_segments "${overlap_segments}" \\
        --bimfile "${bim}" \\
        --chromosome ${chr} \\
        --out_dir . \\
        --out_prefix "${specificity_id}" \\
        --all_genes
    """
}
