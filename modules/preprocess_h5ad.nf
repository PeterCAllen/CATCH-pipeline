// modules/local/preprocess_h5ad.nf

process PREPROCESS_H5AD {
    tag "$h5ad_raw"
    label 'high_mem'
    container "${projectDir}/environments/py-r-cepo-scdrs.sif"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/scanpy:1.9.3--pyhdfd78af_0' :
    //     'biocontainers/scanpy:1.9.3--pyhdfd78af_0' }"

    input:
    path h5ad_raw

    output:
    path "filtered.h5ad", emit: h5ad

    script:
    """
    python3 -u ${projectDir}/bin/subset_h5ad.py \\
        --input_h5ad ${h5ad_raw} \\
        --output_h5ad filtered.h5ad \\
        --cell_type_col ${params.cell_type_column}
    """
}