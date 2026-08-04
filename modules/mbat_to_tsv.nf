// modules/mbat_to_tsv.nf
// Convert mBAT-combo results to the scDRS z-score table.
// The z-score column is named after the trait because `scdrs munge-gs` uses each
// non-GENE column name as the trait name in the resulting .gs file.

process MBAT_TO_TSV {
    tag "${params.gwas_name}"
    label 'medium_mem'
    publishDir "${params.outdir}/scdrs/gwas", mode: 'copy'

    container "${projectDir}/environments/py-r-cepo-scdrs.sif"

    input:
    path mbat_combined

    output:
    path "${params.gwas_name}.tsv", emit: tsv_file
    path "*.log",                   emit: log

    script:
    """
    Rscript - "${mbat_combined}" "${params.gwas_name}.tsv" "${params.scdrs_top_genes}" "${params.gwas_name}" <<'RSCRIPT' 2>&1 | tee ${params.gwas_name}_mbat_to_tsv.log
    suppressPackageStartupMessages(library(data.table))

    args       <- commandArgs(trailingOnly = TRUE)
    mbat_file  <- args[1]
    tsv_file   <- args[2]
    top_n      <- as.integer(args[3])
    trait_name <- args[4]

    dt <- fread(mbat_file)

    for (col in c("Gene", "P_mBATcombo")) {
      if (!col %in% names(dt)) {
        stop("mBAT output missing '", col, "'. Columns: ", paste(names(dt), collapse = ", "))
      }
    }

    dt <- dt[!is.na(P_mBATcombo)]

    min_nonzero_p <- min(dt\$P_mBATcombo[dt\$P_mBATcombo > 0], na.rm = TRUE)
    dt[P_mBATcombo <= 0, P_mBATcombo := min_nonzero_p / 10]
    dt[P_mBATcombo > 1 - 1e-16, P_mBATcombo := 1 - 1e-16]

    dt_sorted <- dt[order(P_mBATcombo)][1:min(.N, top_n)]

    z <- qnorm(dt_sorted\$P_mBATcombo / 2, lower.tail = FALSE) * sign(0.5 - dt_sorted\$P_mBATcombo)

    out <- data.table(GENE = dt_sorted\$Gene, Z = z)
    setnames(out, "Z", trait_name)

    fwrite(out, tsv_file, sep = "\\t")

    cat(sprintf("Wrote %s: %d genes, trait column '%s'\\n", tsv_file, nrow(out), trait_name))
    cat(sprintf("  Z range: %.3f to %.3f (mean %.3f)\\n", min(z), max(z), mean(z)))
    cat(sprintf("  P range: %.3e to %.3e\\n",
                min(dt_sorted\$P_mBATcombo), max(dt_sorted\$P_mBATcombo)))

    if (any(!is.finite(z))) stop("Non-finite z-scores produced.")
    if (median(z) < 0) {
      stop("Median z-score is negative; the p-to-z sign convention is inverted.")
    }
    RSCRIPT
    """
}
