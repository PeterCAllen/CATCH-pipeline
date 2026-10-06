// modules/conldsc/make_cts_file.nf
// CELLECT rule: make_cts_file

process MAKE_CTS_FILE {
    tag "${specificity_id}"
    label 'low_mem'
    publishDir { "${params.outdir}/conldsc/${specificity_id}/precomputation" }, mode: 'copy'

    container "${projectDir}/environments/cellect-py3.sif"

    input:
    tuple val(specificity_id), path(annotations), path(per_annotation_files)

    output:
    tuple val(specificity_id), path("${specificity_id}.ldcts.txt"), emit: ldcts

    script:
    """
    python3 -u "${projectDir}/bin/cellect_make_cts_file.py" \\
        --annotations "${annotations}" \\
        --run_prefix "${specificity_id}" \\
        --ldscore_dir . \\
        --out "${specificity_id}.ldcts.txt"

    while IFS=\$'\\t' read -r name prefix; do
        for chr in \$(seq 1 22); do
            if [ ! -s "\${prefix}\${chr}.l2.ldscore.gz" ]; then
                echo "ERROR: missing \${prefix}\${chr}.l2.ldscore.gz referenced by the .ldcts file"
                exit 1
            fi
        done
    done < "${specificity_id}.ldcts.txt"

    echo "All .ldcts entries resolve to 22 chromosomes of LD scores."
    """
}
