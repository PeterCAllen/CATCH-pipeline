// subworkflows/scdrs.nf
//
// scDRS driven by mBAT-combo gene scores. Consumes the shared mBAT-combo output
// so GCTA runs once for both this branch and seismic.

include { MBAT_TO_TSV          } from '../modules/mbat_to_tsv'
include { MUNGE_SCDRS_GENESET  } from '../modules/munge_scdrs_geneset'
include { RUN_SCDRS_SCORE      } from '../modules/run_scdrs_score'
include { RUN_SCDRS_DOWNSTREAM } from '../modules/run_scdrs_downstream'

workflow SCDRS {
    take:
        mbat_combined   // combined mBAT-combo results (Gene, P_mBATcombo)
        h5ad            // normalized h5ad
        cov_file        // covariate file, or NO_FILE

    main:
        MBAT_TO_TSV(mbat_combined)
        MUNGE_SCDRS_GENESET(MBAT_TO_TSV.out.tsv_file)

        ch_score_input = h5ad
            .combine(MUNGE_SCDRS_GENESET.out.geneset)
            .combine(cov_file)

        RUN_SCDRS_SCORE(ch_score_input)

        ch_downstream_input = h5ad
            .combine(RUN_SCDRS_SCORE.out.score_file)
            .combine(cov_file)

        RUN_SCDRS_DOWNSTREAM(ch_downstream_input)

        ch_group_results = RUN_SCDRS_DOWNSTREAM.out.group_results
            .flatten()
            .filter { it.name.contains('.scdrs_group') }

    emit:
        tsv_file      = MBAT_TO_TSV.out.tsv_file
        geneset       = MUNGE_SCDRS_GENESET.out.geneset
        score_file    = RUN_SCDRS_SCORE.out.score_file
        group_results = ch_group_results
}
