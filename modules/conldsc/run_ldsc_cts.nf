// modules/conldsc/run_ldsc_cts.nf
// CELLECT rule: prioritize_annotations

process RUN_LDSC_CTS {
    tag "${specificity_id}"
    label 'medium_mem'
    publishDir { "${params.outdir}/conldsc/${specificity_id}/prioritization" }, mode: 'copy'

    container "${projectDir}/environments/ldsc-py3.sif"

    input:
    tuple val(specificity_id),
          path(gwas_munged),
          path(ldcts),
          path(per_annotation_files),
          path(all_genes_files),
          path(baseline_dir),
          path(weights_dir)

    output:
    tuple val(specificity_id), path("${specificity_id}__${params.gwas_name}.cell_type_results.txt"), emit: cell_type_results
    path "*.log", emit: log

    script:
    def baseline_with_dot = params.ref_hg19_baseline.substring(params.ref_hg19_baseline.lastIndexOf('/') + 1)
    def baseline_name     = baseline_with_dot.endsWith('.') ? baseline_with_dot[0..-2] : baseline_with_dot
    def weights_with_dot  = params.ref_hg19_weights.substring(params.ref_hg19_weights.lastIndexOf('/') + 1)
    def weights_name      = weights_with_dot.endsWith('.') ? weights_with_dot[0..-2] : weights_with_dot
    def out_prefix        = "${specificity_id}__${params.gwas_name}"
    """
    for chr in \$(seq 1 22); do
        for ext in annot.gz l2.ldscore.gz l2.M l2.M_5_50; do
            src="${baseline_dir}/${baseline_name}.\${chr}.\${ext}"
            if [ -f "\$src" ]; then
                ln -sf "\$src" .
            else
                echo "ERROR: missing baseline file \$src"
                exit 1
            fi
        done

        wsrc="${weights_dir}/${weights_name}.\${chr}.l2.ldscore.gz"
        if [ -f "\$wsrc" ]; then
            ln -sf "\$wsrc" .
        else
            echo "ERROR: missing weights file \$wsrc"
            exit 1
        fi
    done

    ldsc.py \\
        --h2-cts "${gwas_munged}" \\
        --ref-ld-chr "${baseline_name}.,all_genes_in_${specificity_id}." \\
        --w-ld-chr "${weights_name}." \\
        --ref-ld-chr-cts "${ldcts}" \\
        --out "${out_prefix}" \\
        2>&1 | tee "${out_prefix}.h2cts.stdout.log"

    if [ ! -s "${out_prefix}.cell_type_results.txt" ]; then
        echo "ERROR: ldsc.py --h2-cts produced no ${out_prefix}.cell_type_results.txt"
        exit 1
    fi

    N_IN=\$(wc -l < "${ldcts}")
    N_OUT=\$(( \$(wc -l < "${out_prefix}.cell_type_results.txt") - 1 ))
    echo "Cell types submitted: \$N_IN, returned: \$N_OUT"
    if [ "\$N_OUT" -lt "\$N_IN" ]; then
        echo "WARNING: \$(( N_IN - N_OUT )) cell types were dropped by LDSC. Check the log for '*CTS ERROR*'."
        grep -i 'CTS ERROR' "${out_prefix}.log" || true
    fi
    """
}
