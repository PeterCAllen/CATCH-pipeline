// Diagnostic: verify the Part A fix (intersect mBAT genes with the h5ad's gene
// list before ranking/truncating to top-N) against the already-computed mBAT
// output and our own normalized h5ad from the T2D_06_v3 run. Reuses cached
// mBAT results -- no need to rerun GCTA.
nextflow.enable.dsl = 2

include { SCDRS } from './subworkflows/scdrs'

workflow {
    ch_mbat  = Channel.fromPath(params.mbat_combined, checkIfExists: true)
    ch_h5ad  = Channel.fromPath(params.norm_h5ad,     checkIfExists: true)
    ch_cov   = params.scdrs_cov_file
        ? Channel.fromPath(params.scdrs_cov_file, checkIfExists: true)
        : Channel.value(file("NO_FILE"))

    SCDRS(ch_mbat, ch_h5ad, ch_cov)

    SCDRS.out.geneset.view { f -> "geneset: ${f}" }
    SCDRS.out.group_results.view { f -> "group_results: ${f}" }
}
