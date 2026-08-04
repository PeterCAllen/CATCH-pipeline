// subworkflows/mbat.nf
//
// GCTA mBAT-combo gene-based association test (Li et al. Methods: TSS/TES +/- 10 kb,
// 1000G Phase 3 EUR). Run once and shared by the seismic and scDRS branches.

include { FORMAT_GWAS_FOR_MBAT } from '../modules/format_gwas_for_mbat'
include { PREPARE_MBAT_GENES   } from '../modules/prepare_mbat_genes'
include { RUN_MBAT             } from '../modules/run_mbat'
include { COMBINE_MBAT         } from '../modules/combine_mbat'

workflow MBAT {
    take:
        gwas_sumstats   // raw GWAS summary statistics
        gene_coords     // CEPO/scDRS gene coordinates (with header)
        genome_build    // hg19 or hg38

    main:
        def plink_dir    = params.ref_hg19_plink_dir
        def plink_prefix = new File(params.ref_hg19_plink_prefix).name

        ch_ref_bim = Channel.value(file(params.ref_hg19_bim_file, checkIfExists: true))

        FORMAT_GWAS_FOR_MBAT(gwas_sumstats, genome_build, ch_ref_bim)
        PREPARE_MBAT_GENES(gene_coords)

        ch_chr_plink = Channel.of(1..22).map { chr ->
            tuple(chr, plink_prefix, [
                file("${plink_dir}/${plink_prefix}.${chr}.bed", checkIfExists: true),
                file("${plink_dir}/${plink_prefix}.${chr}.bim", checkIfExists: true),
                file("${plink_dir}/${plink_prefix}.${chr}.fam", checkIfExists: true)
            ])
        }

        ch_mbat_input = Channel.of(1..22)
            .combine(FORMAT_GWAS_FOR_MBAT.out.formatted_gwas)
            .combine(PREPARE_MBAT_GENES.out.mbat_genes)
            .combine(ch_chr_plink, by: 0)
            .map { chr, gwas, genes, prefix, plink_files ->
                tuple(chr, gwas, genes, prefix, plink_files)
            }

        RUN_MBAT(ch_mbat_input)

        ch_mbat_collected = RUN_MBAT.out.mbat_result
            .map { chr, f -> f }
            .collect()
            .map { files ->
                if (files.size() != 22) {
                    error "mBAT produced ${files.size()} chromosome files, expected 22."
                }
                def gwas_prefix = files[0].baseName.replaceAll(~/_chr\d+\.gene\.assoc$/, '')
                tuple(gwas_prefix, files)
            }

        COMBINE_MBAT(ch_mbat_collected)

    emit:
        formatted_gwas = FORMAT_GWAS_FOR_MBAT.out.formatted_gwas
        mbat_genes     = PREPARE_MBAT_GENES.out.mbat_genes
        mbat_combined  = COMBINE_MBAT.out.mbat_combined
}
