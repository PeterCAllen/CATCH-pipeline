// subworkflows/mbat.nf
//
// GCTA mBAT-combo gene-based association test (Li et al. Methods: TSS/TES +/- 10 kb,
// 1000G Phase 3 EUR). Run once and shared by the seismic and scDRS branches.
//
// LD reference is MAGMA's g1000_eur (503 EUR, ~22.6M SNPs, one combined whole-genome
// PLINK set, not the HapMap3-intersected panel conLDSC/LDSC uses) -- mBAT-combo's test
// statistic depends on SVD of the local LD matrix within each gene window, so a denser,
// purpose-built reference matters here in a way it doesn't for genome-wide LDSC.

include { FORMAT_GWAS_FOR_MBAT } from '../modules/format_gwas_for_mbat'
include { RUN_MBAT             } from '../modules/run_mbat'
include { COMBINE_MBAT         } from '../modules/combine_mbat'

workflow MBAT {
    take:
        gwas_sumstats   // raw GWAS summary statistics
        mbat_gene_coords // mBAT gene list from PREPARE_GENE_COORDS.out.mbat_loc: chr, start, end,
                         // Gene -- no header, unfiltered, already sorted by chr/start/end
        genome_build    // hg19 or hg38

    main:
        def plink_prefix = new File(params.ref_hg19_mbat_plink_prefix).name

        if (params.gwas_cojo) {
            // Pre-formatted GCTA-COJO .ma supplied directly; skip FORMAT_GWAS_FOR_MBAT.
            ch_formatted_gwas = Channel.fromPath(params.gwas_cojo, checkIfExists: true)
        } else {
            // rsID lookup must use the same panel mBAT-combo will run against below --
            // otherwise SNPs outside the HapMap3-intersected panel never get an rsID and
            // are dropped before GCTA ever sees them, capping mBAT's SNP set to the
            // sparser panel regardless of which --bfile it's pointed at.
            ch_ref_bim = Channel.value(file("${params.ref_hg19_mbat_plink_prefix}.bim", checkIfExists: true))
            FORMAT_GWAS_FOR_MBAT(gwas_sumstats, genome_build, ch_ref_bim)
            ch_formatted_gwas = FORMAT_GWAS_FOR_MBAT.out.formatted_gwas
        }

        // One combined whole-genome bfile, broadcast to all 22 --chr-restricted tasks
        // (g1000_eur ships as a single bed/bim/fam set, not split per chromosome).
        ch_mbat_plink = Channel.value(tuple(plink_prefix, [
            file("${params.ref_hg19_mbat_plink_prefix}.bed", checkIfExists: true),
            file("${params.ref_hg19_mbat_plink_prefix}.bim", checkIfExists: true),
            file("${params.ref_hg19_mbat_plink_prefix}.fam", checkIfExists: true)
        ]))

        ch_mbat_input = Channel.of(1..22)
            .combine(ch_formatted_gwas)
            .combine(mbat_gene_coords)
            .combine(ch_mbat_plink)
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
        formatted_gwas = ch_formatted_gwas
        mbat_genes     = mbat_gene_coords
        mbat_combined  = COMBINE_MBAT.out.mbat_combined
}
