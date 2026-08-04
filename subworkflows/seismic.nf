// subworkflows/seismic.nf
//
// seismic-mBAT-combo. Li et al. Methods: "In the structure of Seismic, the disease genes
// z-statistics could be replaced with mBAT-combo based z-statistics".
// seismic computes its own specificity via calc_specificity(); it does not consume Cepo or GES.

include { MBAT_TO_SEISMIC_ZSTAT } from '../modules/mbat_to_seismic_zstat'
include { RUN_SEISMIC           } from '../modules/run_seismic'

workflow SEISMIC {
    take:
        mbat_combined   // combined mBAT-combo results (Gene, P_mBATcombo)
        h5ad            // normalized h5ad
        gene_coords     // CEPO/scDRS gene coordinates (Gene, gene_name)

    main:
        MBAT_TO_SEISMIC_ZSTAT(mbat_combined, gene_coords)

        ch_seismic_input = h5ad.combine(MBAT_TO_SEISMIC_ZSTAT.out.zstat)

        RUN_SEISMIC(ch_seismic_input)

    emit:
        zstat       = MBAT_TO_SEISMIC_ZSTAT.out.zstat
        results     = RUN_SEISMIC.out.results
        specificity = RUN_SEISMIC.out.specificity
}
