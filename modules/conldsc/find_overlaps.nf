// modules/conldsc/find_overlaps.nf
// CELLECT rule: find_overlaps

process FIND_OVERLAPS {
    tag "chr${chr}"
    label 'medium_mem'
    publishDir "${params.outdir}/conldsc/shared/bed", mode: 'copy'

    container "${projectDir}/environments/cellect-py3.sif"

    input:
    tuple val(chr), path(gene_bed)

    output:
    tuple val(chr), path("overlap_segments_${params.conldsc_window_kb}kb.${chr}.bed"), emit: overlap_segments

    script:
    """
    bedops --partition "${gene_bed}" \\
        | bedmap --echo --echo-map-id-uniq --delim '\\t' - "${gene_bed}" \\
        > "overlap_segments_${params.conldsc_window_kb}kb.${chr}.bed"

    if [ ! -s "overlap_segments_${params.conldsc_window_kb}kb.${chr}.bed" ]; then
        echo "ERROR: empty overlap segments for chr${chr}"
        exit 1
    fi
    """
}
