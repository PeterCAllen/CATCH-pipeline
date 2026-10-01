// modules/mbat_to_tsv.nf
// Convert mBAT-combo results to the scDRS z-score table.
// The z-score column is named after the trait because `scdrs munge-gs` uses each
// non-GENE column name as the trait name in the resulting .gs file.
//
// Intersects mBAT genes with the h5ad's actual gene list BEFORE ranking/truncating
// to the top N -- a gene selected as "top 1000" that isn't present in the h5ad
// scDRS will score against is useless and silently shrinks the final gene set.

process MBAT_TO_TSV {
    tag "${params.gwas_name}"
    label 'medium_mem'
    publishDir "${params.outdir}/scdrs/gwas", mode: 'copy'

    container "${projectDir}/environments/py-r-cepo-scdrs.sif"

    input:
    path mbat_combined
    path h5ad

    output:
    path "${params.gwas_name}.tsv", emit: tsv_file
    path "*.log",                   emit: log

    script:
    """
    python3 -c "
    import anndata
    import pandas as pd
    ad = anndata.read_h5ad('${h5ad}', backed='r')
    pd.Series(ad.var_names, name='Gene').to_csv('gene_list.csv', index=False)
    "

    Rscript - "${mbat_combined}" "gene_list.csv" "${params.gwas_name}.tsv" "${params.scdrs_top_genes}" "${params.gwas_name}" <<'RSCRIPT' 2>&1 | tee ${params.gwas_name}_mbat_to_tsv.log
    suppressPackageStartupMessages(library(data.table))

    args           <- commandArgs(trailingOnly = TRUE)
    mbat_file      <- args[1]
    gene_list_file <- args[2]
    tsv_file       <- args[3]
    top_n          <- as.integer(args[4])
    trait_name     <- args[5]

    dt <- fread(mbat_file)

    for (col in c("Gene", "P_mBATcombo")) {
      if (!col %in% names(dt)) {
        stop("mBAT output missing '", col, "'. Columns: ", paste(names(dt), collapse = ", "))
      }
    }

    dt <- dt[!is.na(P_mBATcombo)]

    gene_list <- fread(gene_list_file)
    n_before  <- nrow(dt)
    dt        <- dt[Gene %in% gene_list\$Gene]
    cat(sprintf("Gene-list intersection: %d -> %d mBAT genes present in the h5ad\\n", n_before, nrow(dt)))

    min_nonzero_p <- min(dt\$P_mBATcombo[dt\$P_mBATcombo > 0], na.rm = TRUE)
    dt[P_mBATcombo <= 0, P_mBATcombo := min_nonzero_p]
    dt[P_mBATcombo > 1 - 1e-16, P_mBATcombo := 1 - 1e-16]

    dt_sorted <- dt[order(P_mBATcombo)][1:min(.N, top_n)]

    z <- qnorm(pmax(dt_sorted\$P_mBATcombo, .Machine\$double.xmin) / 2, lower.tail = FALSE)

    out <- data.table(GENE = dt_sorted\$Gene, Z = z)
    setnames(out, "Z", trait_name)

    fwrite(out, tsv_file, sep = "\\t")

    cat(sprintf("Wrote %s: %d genes, trait column '%s'\\n", tsv_file, nrow(out), trait_name))
    cat(sprintf("  Z range: %.3f to %.3f (mean %.3f)\\n", min(z), max(z), mean(z)))
    cat(sprintf("  P range: %.3e to %.3e\\n",
                min(dt_sorted\$P_mBATcombo), max(dt_sorted\$P_mBATcombo)))

    if (any(!is.finite(z))) stop("Non-finite z-scores produced.")
    if (any(z < 0)) stop("Negative z-score produced; p-to-z conversion is broken.")
    RSCRIPT
    """
}
