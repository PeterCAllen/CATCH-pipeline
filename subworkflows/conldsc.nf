// subworkflows/conldsc.nf
//
// Native Nextflow port of CELLECT-LDSC (perslab/CELLECT, cellect-ldsc.snakefile).
// Runs once per specificity matrix (e.g. cepo_norm, ges) and emits CELLECT's
// prioritization.csv schema: gwas, specificity_id, annotation, beta, beta_se, pvalue.

include { MAKE_ALL_GENES_BACKGROUND    } from '../modules/conldsc/make_all_genes'
include { FORMAT_GENES                 } from '../modules/conldsc/format_genes'
include { FIND_OVERLAPS                } from '../modules/conldsc/find_overlaps'
include { MAKE_ANNOT                   } from '../modules/conldsc/make_annot'
include { MAKE_ANNOT_ALL_GENES         } from '../modules/conldsc/make_annot'
include { COMPUTE_LD_SCORES            } from '../modules/conldsc/compute_ld_scores'
include { COMPUTE_LD_SCORES_ALL_GENES  } from '../modules/conldsc/compute_ld_scores'
include { SPLIT_LD_SCORES              } from '../modules/conldsc/split_ld_scores'
include { MAKE_CTS_FILE                } from '../modules/conldsc/make_cts_file'
include { RUN_LDSC_CTS                 } from '../modules/conldsc/run_ldsc_cts'
include { PARSE_LDSC_RESULTS           } from '../modules/conldsc/parse_results'

workflow CONLDSC {
    take:
        specificity     // tuple(specificity_id, sanitized specificity csv)
        annotations     // tuple(specificity_id, annotation name list)
        gwas_munged     // LDSC .sumstats.gz
        cellect_coords  // CELLECT-format gene coordinates (6 cols, no header)

    main:
        def plink_dir    = params.ref_hg19_plink_dir
        def plink_prefix = new File(params.ref_hg19_plink_prefix).name
        def baseline_dir = params.ref_hg19_baseline.substring(0, params.ref_hg19_baseline.lastIndexOf('/'))
        def weights_dir  = params.ref_hg19_weights.substring(0, params.ref_hg19_weights.lastIndexOf('/'))

        ch_chr_sizes  = Channel.value(file(params.ref_hg19_chr_sizes,  checkIfExists: true))
        ch_print_snps = Channel.value(file(params.ref_hg19_print_snps, checkIfExists: true))
        ch_baseline   = Channel.value(file(baseline_dir, checkIfExists: true))
        ch_weights    = Channel.value(file(weights_dir,  checkIfExists: true))

        ch_chr_bim = Channel.of(1..22).map { chr ->
            tuple(chr, file("${plink_dir}/${plink_prefix}.${chr}.bim", checkIfExists: true))
        }

        ch_chr_plink = Channel.of(1..22).map { chr ->
            tuple(chr, [
                file("${plink_dir}/${plink_prefix}.${chr}.bed", checkIfExists: true),
                file("${plink_dir}/${plink_prefix}.${chr}.bim", checkIfExists: true),
                file("${plink_dir}/${plink_prefix}.${chr}.fam", checkIfExists: true)
            ])
        }

        // --- Shared across specificity matrices: windowed gene BEDs and overlap segments ---

        FORMAT_GENES(cellect_coords, ch_chr_sizes)

        ch_gene_beds = FORMAT_GENES.out.gene_beds
            .flatten()
            .map { bed ->
                def m = (bed.name =~ /\.(\d{1,2})\.bed$/)
                if (!m) error "Cannot parse chromosome from gene BED filename: ${bed.name}"
                tuple(m[0][1] as Integer, bed)
            }

        FIND_OVERLAPS(ch_gene_beds)

        ch_segments_bim = FIND_OVERLAPS.out.overlap_segments
            .join(ch_chr_bim, by: 0)      // tuple(chr, segments, bim)

        // --- Per specificity matrix ---

        MAKE_ALL_GENES_BACKGROUND(specificity)

        ch_annot_input = specificity
            .combine(ch_segments_bim)
            .map { spec_id, matrix, chr, segments, bim ->
                tuple(spec_id, chr, matrix, segments, bim)
            }

        MAKE_ANNOT(ch_annot_input)

        ch_all_genes_input = MAKE_ALL_GENES_BACKGROUND.out.all_genes
            .combine(ch_segments_bim)
            .map { spec_id, all_genes, chr, segments, bim ->
                tuple(spec_id, chr, all_genes, segments, bim)
            }

        MAKE_ANNOT_ALL_GENES(ch_all_genes_input)

        // --- LD scores ---

        ch_ldsc_input = MAKE_ANNOT.out.annot
            .map { spec_id, chr, annot -> tuple(chr, spec_id, annot) }
            .combine(ch_chr_plink, by: 0)
            .combine(ch_print_snps)
            .map { chr, spec_id, annot, plink, print_snps ->
                tuple(spec_id, chr, annot, plink, print_snps)
            }

        COMPUTE_LD_SCORES(ch_ldsc_input)

        ch_ldsc_all_genes_input = MAKE_ANNOT_ALL_GENES.out.annot
            .map { spec_id, chr, annot -> tuple(chr, spec_id, annot) }
            .combine(ch_chr_plink, by: 0)
            .combine(ch_print_snps)
            .map { chr, spec_id, annot, plink, print_snps ->
                tuple(spec_id, chr, annot, plink, print_snps)
            }

        COMPUTE_LD_SCORES_ALL_GENES(ch_ldsc_all_genes_input)

        // --- Split the combined LD scores into one set per cell type ---

        SPLIT_LD_SCORES(COMPUTE_LD_SCORES.out.ldscores)

        ch_per_annotation = SPLIT_LD_SCORES.out.per_annotation
            .map { spec_id, files -> tuple(spec_id, files) }
            .groupTuple(by: 0)
            .map { spec_id, file_lists -> tuple(spec_id, file_lists.flatten()) }

        ch_cts_input = annotations.join(ch_per_annotation, by: 0)

        MAKE_CTS_FILE(ch_cts_input)

        // --- Regression: one ldsc.py --h2-cts job per specificity matrix ---

        ch_all_genes_grouped = COMPUTE_LD_SCORES_ALL_GENES.out.ldscores
            .map { spec_id, ldscore, m, m550, annot -> tuple(spec_id, [ldscore, m, m550, annot]) }
            .groupTuple(by: 0)
            .map { spec_id, file_lists -> tuple(spec_id, file_lists.flatten()) }

        ch_cts_run_input = MAKE_CTS_FILE.out.ldcts
            .join(ch_per_annotation,   by: 0)
            .join(ch_all_genes_grouped, by: 0)
            .combine(gwas_munged)
            .combine(ch_baseline)
            .combine(ch_weights)
            .map { spec_id, ldcts, per_annot, all_genes, munged, baseline, weights ->
                tuple(spec_id, munged, ldcts, per_annot, all_genes, baseline, weights)
            }

        RUN_LDSC_CTS(ch_cts_run_input)

        PARSE_LDSC_RESULTS(RUN_LDSC_CTS.out.cell_type_results)

    emit:
        prioritization    = PARSE_LDSC_RESULTS.out.prioritization   // tuple(specificity_id, csv)
        cell_type_results = RUN_LDSC_CTS.out.cell_type_results
}
