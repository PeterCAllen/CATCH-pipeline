// Standalone diagnostic: run scDRS directly on the h5ad's native .X (Seurat
// SCTransform-normalized, 17199 genes) instead of our own raw-promoted +
// CPM-log2-renormalized matrix. Ang's scdrs job script points scDRS straight
// at "${DATASET_NAME}.h5ad" (the same file she sent us) with
// --flag-raw-count False, which is consistent with her having fed scDRS this
// native X directly rather than a re-normalized version. Reuses the mBAT-derived
// gene set (.gs) from the T2D_06_v3 run since that only depends on mBAT, not
// on h5ad normalization.
nextflow.enable.dsl = 2

include { RUN_SCDRS_SCORE      } from './modules/run_scdrs_score'
include { RUN_SCDRS_DOWNSTREAM } from './modules/run_scdrs_downstream'

workflow {
    ch_h5ad    = Channel.fromPath(params.native_h5ad,  checkIfExists: true)
    ch_geneset = Channel.fromPath(params.scdrs_geneset, checkIfExists: true)
    ch_cov     = params.scdrs_cov_file
        ? Channel.fromPath(params.scdrs_cov_file, checkIfExists: true)
        : Channel.value(file("NO_FILE"))

    ch_score_input = ch_h5ad
        .combine(ch_geneset)
        .combine(ch_cov)

    RUN_SCDRS_SCORE(ch_score_input)

    ch_downstream_input = ch_h5ad
        .combine(RUN_SCDRS_SCORE.out.score_file)
        .combine(ch_cov)

    RUN_SCDRS_DOWNSTREAM(ch_downstream_input)

    RUN_SCDRS_DOWNSTREAM.out.group_results.view { f -> "group_results: ${f}" }
}
