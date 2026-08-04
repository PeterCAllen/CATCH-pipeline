// subworkflows/metrics.nf
//
// Cell-type specificity metrics for the CATCH conLDSC branches:
//   Cepo   -> cepo_norm (raw stats) and cepo_s (esw_star rank-normalized)
//   CELLEX -> ges
//
// Emits a channel of tuple(specificity_id, sanitized CELLECT-ready matrix).

include { asBool } from '../lib/util'

include { RUN_CEPO             } from '../modules/run_cepo'
include { CEPO_ESW_STAR        } from '../modules/cepo_esw_star'
include { RUN_CELLEX           } from '../modules/run_cellex'
include { SANITIZE_SPECIFICITY } from '../modules/sanitize_specificity'

workflow METRICS {
    take:
        h5ad        // normalized h5ad (log2(TPM+1), sanitized cell type labels)

    main:
        ch_raw          = Channel.empty()
        ch_cepo_rds     = Channel.empty()
        ch_cepo_pvalues = Channel.empty()
        ch_ges          = Channel.empty()

        if (asBool(params.run_cepo)) {
            if (!(params.conldsc_cepo_variant in ['cepo_norm', 'cepo_s'])) {
                error "Unknown --conldsc_cepo_variant '${params.conldsc_cepo_variant}'. Use 'cepo_norm' or 'cepo_s'."
            }

            RUN_CEPO(h5ad)
            ch_cepo_rds     = RUN_CEPO.out.rds
            ch_cepo_pvalues = RUN_CEPO.out.pvalues

            if (params.conldsc_cepo_variant == 'cepo_s') {
                CEPO_ESW_STAR(RUN_CEPO.out.stats, RUN_CEPO.out.pvalues)
                ch_raw = ch_raw.mix(CEPO_ESW_STAR.out.stats.map { m -> tuple('cepo_s', m, 'comma') })
            }
            else {
                ch_raw = ch_raw.mix(RUN_CEPO.out.stats.map { m -> tuple('cepo_norm', m, 'comma') })
            }
        }

        if (asBool(params.run_cellex)) {
            RUN_CELLEX(h5ad)
            ch_ges = RUN_CELLEX.out.ges
            ch_raw = ch_raw.mix(ch_ges.map { m -> tuple('ges', m, 'comma') })
        }

        SANITIZE_SPECIFICITY(ch_raw)

    emit:
        specificity  = SANITIZE_SPECIFICITY.out.matrix       // tuple(specificity_id, csv)
        annotations  = SANITIZE_SPECIFICITY.out.annotations  // tuple(specificity_id, txt)
        cepo_rds     = ch_cepo_rds
        cepo_pvalues = ch_cepo_pvalues
        ges          = ch_ges
}
