// modules/local/make_ldsc_annots.nf

process MAKE_LDSC_ANNOTS {
    tag "$bed"
    label 'low_mem'
    conda (params.enable_conda ? "${projectDir}/environments/ldsctools.yml" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ldsc:1.0.1--py27h7814781_4' :
        'biocontainers/ldsc:1.0.1--py27h7814781_4' }"

    input:
    tuple path(bed), path(bim)
    path ldsc_ref_dir

    output:
    path("*.annot.gz"), emit: annot_file
    path("*.l2.ldscore.gz"), emit: ldscore_file
    tuple val(celltype), val(chr), path("*.l2.M"), emit: l2_m
    tuple val(celltype), val(chr), path("*.l2.M_5_50"), emit: l2_m_5_50

    script:
    def chr = bed.name.replaceAll(~/.+\.(\d+)\.bed$/, '$1')
    def celltype = bed.baseName.replaceAll(~/\.\d+$/, '')
    """
    python "${projectDir}/bin/make_annot.py" \\
      --bed-file "${bed}" \\
      --bimfile "${bim}" \\
      --annot-file "${celltype}.${chr}.annot.gz"

    ldsc.py \\
        --l2 \\
        --bfile "${ldsc_ref_dir}/plink_files/1000G.EUR.hg38.${chr}" \\
        --ld-wind-cm 1 \\
        --thin-annot \\
        --annot "${celltype}.${chr}.annot.gz" \\
        --out "${celltype}.${chr}" \\
        --print-snps "${params.ref_hg19_hapmap3}"
    """
}
