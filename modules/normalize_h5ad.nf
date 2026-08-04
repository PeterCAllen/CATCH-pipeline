// modules/normalize_h5ad.nf

include { asBool } from '../lib/util'

process NORMALIZE_H5AD {
    tag "${h5ad.simpleName}"
    label 'high_mem'
    publishDir "${params.outdir}/preprocessed", mode: 'copy'

    container "${projectDir}/environments/py-r-cepo-scdrs.sif"

    input:
    path h5ad
    path gene_coords

    output:
    path "${h5ad.simpleName}.catch.h5ad",           emit: h5ad
    path "${h5ad.simpleName}.celltype_mapping.tsv", emit: cell_type_mapping
    path "${h5ad.simpleName}.celltype_counts.tsv",  emit: cell_type_counts
    path "*.log",                                   emit: log

    script:
    def use_raw = asBool(params.h5ad_use_raw) ? '--use_raw' : ''
    """
    python3 -u "${projectDir}/bin/normalize_h5ad.py" \\
        --input_h5ad "${h5ad}" \\
        --output_h5ad "${h5ad.simpleName}.catch.h5ad" \\
        --output_mapping "${h5ad.simpleName}.celltype_mapping.tsv" \\
        --output_cell_counts "${h5ad.simpleName}.celltype_counts.tsv" \\
        --cell_type_col "${params.cell_type_col}" \\
        --gene_coords "${gene_coords}" \\
        --min_cells ${params.min_cells_per_celltype} \\
        --input_scale ${params.h5ad_input_scale} \\
        --log_base ${params.h5ad_log_base} \\
        ${use_raw} \\
        2>&1 | tee ${h5ad.simpleName}.normalize_h5ad.log
    """
}
