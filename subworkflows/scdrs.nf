// subworkflows/scdrs.nf

//
// Subworkflow for scDRS Analysis (Single-cell Disease Relevance Score)
// Genome-build aware, follows the complete mBAT + scDRS workflow
//

// Import the modules required for this subworkflow
include { FORMAT_GWAS_FOR_SCDRS   } from '../modules/format_gwas_for_scdrs'
include { PREPARE_MBAT_GENES      } from '../modules/prepare_mbat_genes'
include { RUN_MBAT                } from '../modules/run_mbat'
include { COMBINE_MBAT            } from '../modules/combine_mbat'
include { MBAT_TO_TSV             } from '../modules/mbat_to_tsv'
include { MUNGE_SCDRS_GENESET     } from '../modules/munge_scdrs_geneset'
include { RUN_SCDRS_SCORE         } from '../modules/run_scdrs_score'
include { RUN_SCDRS_DOWNSTREAM    } from '../modules/run_scdrs_downstream'

workflow SCDRS {
    take:
        gwas_sumstats       // GWAS summary statistics file (munged or raw with rsIDs)
        h5ad_file           // Single-cell h5ad file
        gene_coords         // Gene coordinates (for CEPO/scDRS)
        genome_build        // hg19 or hg38

    main:
        // Create a channel that selects the appropriate BIM file based on genome build
        ch_ref_bim = genome_build
            .map { build ->
                def ref_bim_path = (build in ['hg38', 'GRCh38']) 
                    ? params.ref_hg38_bim_file
                    : params.ref_hg19_bim_file
                file(ref_bim_path, checkIfExists: true)
            }
        
        // ================================================================================
        // STEP 0: Format GWAS Summary Statistics for mBAT
        // ================================================================================
        // Format: SNP, A1, A2, freq, b, se, p, N
        
        FORMAT_GWAS_FOR_SCDRS(
            gwas_sumstats,
            genome_build,
            ch_ref_bim
        )

        // ================================================================================
        // STEP 0b: Prepare Gene List for mBAT
        // ================================================================================
        // Convert gene coordinates to mBAT format (Gene, Chr, Start, End)
        
        PREPARE_MBAT_GENES(
            gene_coords
        )

        // ================================================================================
        // STEP 1: Run GCTA mBAT per Chromosome
        // ================================================================================
        // Gene-based association test using multivariate Bayesian approach
        
        // Get PLINK reference files based on genome build
        ch_plink_prefix = genome_build.map { build ->
            (build in ['hg38', 'GRCh38']) 
                ? params.ref_hg38_plink_prefix 
                : params.ref_hg19_plink_prefix
        }
        
        // Create channel for each chromosome's PLINK files
        ch_chr_plink = ch_plink_prefix.flatMap { prefix ->
            (1..22).collect { chr ->
                tuple(
                    chr,
                    [
                        file("${prefix}.${chr}.bed", checkIfExists: true),
                        file("${prefix}.${chr}.bim", checkIfExists: true),
                        file("${prefix}.${chr}.fam", checkIfExists: true)
                    ]
                )
            }
        }
        
        // Combine inputs for mBAT: chr, formatted_gwas, mbat_genes, plink_files
        ch_mbat_input = Channel.of(1..22)
            .combine(FORMAT_GWAS_FOR_SCDRS.out.formatted_gwas)
            .combine(PREPARE_MBAT_GENES.out.mbat_genes)
            .combine(ch_chr_plink, by: 0)
            .map { chr, gwas, genes, plink_files ->
                tuple(chr, gwas, genes, plink_files)
            }
        
        RUN_MBAT(ch_mbat_input)

        // ================================================================================
        // STEP 2: Combine mBAT Results from All Chromosomes
        // ================================================================================
        
        // Group mBAT results by GWAS prefix
        // ch_mbat_grouped = RUN_MBAT.out.mbat_result
        //     .map { chr, file ->
        //         def gwas_prefix = file.simpleName.replaceAll(~/_chr\\d+\\.gene\\.assoc$/, '')
        //         tuple(gwas_prefix, file)
        //     }
        //     .groupTuple()
        
        // COMBINE_MBAT(ch_mbat_grouped)

        // ================================================================================
        // STEP 3: Convert mBAT Results to TSV with Z-scores
        // ================================================================================
        // Use ACATO (Cauchy combination) for numerical stability with extreme p-values
        
        // MBAT_TO_TSV(COMBINE_MBAT.out.mbat_combined)

        // ================================================================================
        // STEP 4: Munge to scDRS Gene Set Format
        // ================================================================================
        // Convert TSV with Z-scores to scDRS .gs format
        
        // MUNGE_SCDRS_GENESET(MBAT_TO_TSV.out.tsv_file)

        // ================================================================================
        // STEP 5: Run scDRS Compute-Score
        // ================================================================================
        // Calculate disease relevance scores for each cell
        
        // Combine h5ad with geneset
        // ch_score_input = h5ad_file
        //     .combine(MUNGE_SCDRS_GENESET.out.geneset)
        
        // RUN_SCDRS_SCORE(ch_score_input)

        // ================================================================================
        // STEP 6: Run scDRS Downstream Analysis
        // ================================================================================
        // Cell type association testing
        
        // Combine h5ad with score file
        // ch_downstream_input = h5ad_file
        //     .combine(RUN_SCDRS_SCORE.out.score_file)
        
        // RUN_SCDRS_DOWNSTREAM(ch_downstream_input)

    emit:
        formatted_gwas = FORMAT_GWAS_FOR_SCDRS.out.formatted_gwas
        // mbat_combined  = COMBINE_MBAT.out.mbat_combined
        // tsv_file       = MBAT_TO_TSV.out.tsv_file
        // geneset        = MUNGE_SCDRS_GENESET.out.geneset
        // score_file     = RUN_SCDRS_SCORE.out.score_file
        // group_results  = RUN_SCDRS_DOWNSTREAM.out.group_results
}




