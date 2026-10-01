// Standalone diagnostic: feed Ang's own Cepo output directly into our conLDSC
// branch (bypassing our RUN_CEPO computation entirely) to isolate whether the
// remaining conLDSC-Cepo p-value gap vs her published table is explained by the
// upstream Cepo/h5ad-normalization difference, or whether something downstream
// (conLDSC/LDSC) also diverges.
nextflow.enable.dsl = 2

include { SANITIZE_SPECIFICITY } from './modules/sanitize_specificity'
include { CONLDSC              } from './subworkflows/conldsc'

workflow {
    ch_raw = Channel.of(
        tuple('cepo_norm', file(params.ang_cepo_csv, checkIfExists: true), 'comma')
    )

    SANITIZE_SPECIFICITY(ch_raw)

    ch_gwas_munged    = Channel.fromPath(params.gwas_sumstats_munged, checkIfExists: true)
    ch_cellect_coords = file(params.cellect_coords, checkIfExists: true)

    CONLDSC(
        SANITIZE_SPECIFICITY.out.matrix,
        SANITIZE_SPECIFICITY.out.annotations,
        ch_gwas_munged,
        ch_cellect_coords
    )

    CONLDSC.out.cell_type_results.view { spec_id, f -> "cell_type_results: ${spec_id} -> ${f}" }
    CONLDSC.out.prioritization.view    { spec_id, f -> "prioritization: ${spec_id} -> ${f}" }
}
