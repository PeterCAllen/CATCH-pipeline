// subworkflows/ldsc.nf

//
// Subworkflow for LDSC Partitioned Heritability Analysis (with Continuous Annotations)
// Genome-build aware with proper failsafes
//

// Import the modules required for this subworkflow
include { MUNGE_GWAS_FOR_LDSC         } from '../modules/munge_gwas_for_ldsc'
include { CREATE_CONTINUOUS_BEDS      } from '../modules/create_continuous_beds'
include { CREATE_LDSC_ANNOT           } from '../modules/create_ldsc_annot'
include { COMPUTE_LDSC_SCORES         } from '../modules/compute_ldsc_scores'
include { RUN_SLDSC                   } from '../modules/run_sldsc'
include { COMPUTE_QUANTILE_M          } from '../modules/quantile_analysis'
include { COMPUTE_QUANTILE_H2G        } from '../modules/quantile_analysis'
include { COMPARE_ANNOTATIONS        } from '../modules/quantile_analysis'

workflow LDSC {
    take:
        gwas_sumstats       // GWAS summary statistics file
        cepo_stats          // CEPO results (gene x cell-type matrix)
        ldsc_gene_coords    // Gene coordinates for LDSC (with window)
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
        // STEP 1: Munge GWAS Summary Statistics for LDSC
        // ================================================================================
        // Format: Adds rsIDs, prepares for munge_sumstats.py
        
        MUNGE_GWAS_FOR_LDSC(
            gwas_sumstats,
            genome_build,
            ch_ref_bim
        )

        // ================================================================================
        // STEP 2: Convert CEPO Results to BED Files (one per cell type)
        // ================================================================================
        // Each BED file contains genomic regions with continuous annotation values
        
        CREATE_CONTINUOUS_BEDS(
            cepo_stats,
            ldsc_gene_coords,
            params.ldsc_window_kb
        )

        // ================================================================================
        // STEP 3: Create LDSC Annotations (parallelized by cell type x chromosome)
        // ================================================================================
        // Split into two sub-steps: 3a (R script) and 3b (LDSC Python)
        
        // Extract individual BED files and flatten into separate emissions
        ch_beds = CREATE_CONTINUOUS_BEDS.out.bed_dir
            .flatten()
            .map { bed_file ->
                def cell_type = bed_file.simpleName
                tuple(cell_type, bed_file)
            }
        
        // Determine reference paths based on genome build
        ch_ref_params = genome_build.map { build ->
            def is_hg38 = (build in ['hg38', 'GRCh38'])
            // derive baseline/weights directories from the known ref_base_dir layout
            def baseline_dir = is_hg38 ? "${params.ref_base_dir}/hg38/baseline" : "${params.ref_base_dir}/hg19/baseline"
            def weights_dir  = is_hg38 ? "${params.ref_base_dir}/hg38/weights"  : "${params.ref_base_dir}/hg19/weights"
            [
                plink_dir: is_hg38 ? params.ref_hg38_plink_dir : params.ref_hg19_plink_dir,
                plink_prefix: is_hg38 ? "1000G.EUR.hg38" : "1000G.EUR.hg19",
                hapmap3: is_hg38 ? params.ref_hg38_hapmap3 : params.ref_hg19_hapmap3,
                baseline_dir: baseline_dir,
                weights_dir: weights_dir
            ]
        }
        
        // STEP 3a: Create annotation files (R script in py-r-cepo-scdrs container)
        // Create channels for each chromosome's BIM file only
        ch_chr_bim = ch_ref_params.flatMap { params_map ->
            (1..22).collect { chr ->
                tuple(
                    chr,
                    file("${params_map.plink_dir}/${params_map.plink_prefix}.${chr}.bim", checkIfExists: true)
                )
            }
        }
        
        // Combine BED files with BIM files for annotation creation
        ch_annot_input = ch_beds
            .combine(Channel.of(1..22))
            .map { cell_type, bed, chr ->
                tuple(chr, cell_type, bed)
            }
            .combine(ch_chr_bim, by: 0)
            .map { chr, cell_type, bed, bim ->
                tuple(cell_type, chr, bed, bim)
            }
        
        CREATE_LDSC_ANNOT(
            ch_annot_input.map { cell_type, chr, bed, bim -> tuple(cell_type, chr, bed) },
            genome_build,
            ch_annot_input.map { cell_type, chr, bed, bim -> bim }
        )
        
        // STEP 3b: Compute LD scores (LDSC Python in ldsc_v1.0.1 container)
        // Create channels for each chromosome's full PLINK files (.bed, .bim, .fam)
        ch_chr_plink = ch_ref_params.flatMap { params_map ->
            (1..22).collect { chr ->
                tuple(
                    chr,
                    [
                        file("${params_map.plink_dir}/${params_map.plink_prefix}.${chr}.bed", checkIfExists: true),
                        file("${params_map.plink_dir}/${params_map.plink_prefix}.${chr}.bim", checkIfExists: true),
                        file("${params_map.plink_dir}/${params_map.plink_prefix}.${chr}.fam", checkIfExists: true)
                    ],
                    file(params_map.hapmap3, checkIfExists: true)
                )
            }
        }
        
        // Combine annotations with PLINK files for LD score computation
        ch_ldscore_input = CREATE_LDSC_ANNOT.out.annot_file
            .map { cell_type, chr, annot ->
                tuple(chr, cell_type, annot)
            }
            .combine(ch_chr_plink, by: 0)
            .map { chr, cell_type, annot, plink_files, hapmap3 ->
                tuple(
                    tuple(cell_type, chr, annot),
                    plink_files,
                    hapmap3
                )
            }
        
        COMPUTE_LDSC_SCORES(
            ch_ldscore_input.map { it[0] },  // tuple(cell_type, chr, annot)
            genome_build,
            ch_ldscore_input.map { it[1] },  // plink_files
            ch_ldscore_input.map { it[2] }   // hapmap3
        )
        
        // ================================================================================
        // STEP 3.5: Collect and Verify All LDSC Output Files
        // ================================================================================
        // Group all output files by cell type and ensure all 22 chromosomes succeeded
        
        // Collect annot files from CREATE_LDSC_ANNOT
        ch_annot_grouped = CREATE_LDSC_ANNOT.out.annot_file
            .groupTuple(by: 0)
            .map { cell_type, chr_list, files ->
                def missing = (1..22).findAll { it !in chr_list }
                if (missing) log.warn "Cell type '${cell_type}': Missing annot chromosomes ${missing}"
                tuple(cell_type, files.flatten())
            }
        
        // Collect ldscore files
        ch_ldscore_grouped = COMPUTE_LDSC_SCORES.out.ldscore_file
            .groupTuple(by: 0)
            .map { cell_type, chr_list, files ->
                def missing = (1..22).findAll { it !in chr_list }
                if (missing) log.warn "Cell type '${cell_type}': Missing ldscore chromosomes ${missing}"
                tuple(cell_type, files.flatten())
            }

        // Collect l2.M files
        ch_l2m_grouped = COMPUTE_LDSC_SCORES.out.l2_m
            .groupTuple(by: 0)
            .map { cell_type, chr_list, files ->
                tuple(cell_type, files.flatten())
            }

        // Collect l2.M_5_50 files
        ch_l2m_5_50_grouped = COMPUTE_LDSC_SCORES.out.l2_m_5_50
            .groupTuple(by: 0)
            .map { cell_type, chr_list, files ->
                tuple(cell_type, files.flatten())
            }

        // Combine all annotation files by cell type
        ch_all_annot_files = ch_annot_grouped
            .join(ch_ldscore_grouped, by: 0)
            .join(ch_l2m_grouped, by: 0)
            .join(ch_l2m_5_50_grouped, by: 0)
            .map { cell_type, annot_files, ldscore_files, l2m_files, l2m_5_50_files ->
                tuple(cell_type, annot_files + ldscore_files + l2m_files + l2m_5_50_files)
            }

        // ================================================================================
        // STEP 4: Run Stratified LD Score Regression (one job per cell type)
        // ================================================================================
        
        // Get reference files based on genome build
        ch_ref_files = ch_ref_params.map { params_map ->
            [
                baseline_dir: file(params_map.baseline_dir, checkIfExists: true),
                weights_dir: file(params_map.weights_dir, checkIfExists: true),
                plink_dir: file(params_map.plink_dir, checkIfExists: true),
                plink_prefix: params_map.plink_prefix,
                hapmap3: file(params_map.hapmap3, checkIfExists: true)
            ]
        }

        // Combine all inputs for RUN_SLDSC
        ldsc_input_ch = ch_all_annot_files
            .combine(MUNGE_GWAS_FOR_LDSC.out.gwas_munged)
            .combine(genome_build)
            .combine(ch_ref_files)
            .map { cell_type, annot_files, gwas_munged, build, ref_files ->
                tuple(
                    cell_type,
                    gwas_munged,
                    build,
                    ref_files.baseline_dir,
                    ref_files.weights_dir,
                    ref_files.plink_dir,
                    ref_files.plink_prefix,
                    ref_files.hapmap3,
                    annot_files
                )
            }
        
        RUN_SLDSC(ldsc_input_ch)

        // ================================================================================
        // STEP 5: Quantile-based M Analysis
        // ================================================================================
        
        // STEP 5a: Compute Quantile M values (Perl)
        quantile_m_input_ch = ch_annot_grouped
            .combine(genome_build)
            .combine(ch_ref_files)
            .map { cell_type, annot_files, build, ref_files ->
                tuple(
                    cell_type,
                    annot_files,
                    build,
                    ref_files.baseline_dir,
                    ref_files.plink_dir,
                    ref_files.plink_prefix
                )
            }
        
        COMPUTE_QUANTILE_M(quantile_m_input_ch)

        // STEP 5b: Compute Quantile h2g (R)
        // Create a proper channel from RUN_SLDSC results
        ch_sldsc_results = RUN_SLDSC.out.results
            .map { file -> 
                def cell_type = file.baseName.replaceAll("_sldsc", '')
                tuple(cell_type, file)
            }
        ch_sldsc_parts_delete = RUN_SLDSC.out.delete_vals
            .map { file -> 
                def cell_type = file.baseName.replaceAll("_sldsc", '')
                tuple(cell_type, file)
            }
        
        // Join quantile M output with S-LDSC results by cell type
        quantile_h2g_input_ch = COMPUTE_QUANTILE_M.out.quantile_m
            .join(ch_sldsc_results, by: 0)
            .join(ch_sldsc_parts_delete, by: 0)
        
        COMPUTE_QUANTILE_H2G(quantile_h2g_input_ch)

        // ================================================================================
        // STEP 6: Compare Annotations Across Cell Types
        // ================================================================================
        
        // Collect all quantile results and run comparison
        COMPARE_ANNOTATIONS(
            COMPUTE_QUANTILE_H2G.out.quantile_results.collect(),
            params.gwas_name
        )

    emit:
        munged_gwas = MUNGE_GWAS_FOR_LDSC.out.gwas_munged
        bed_dir     = CREATE_CONTINUOUS_BEDS.out.bed_dir
        results     = RUN_SLDSC.out.results
        quantile_results = COMPUTE_QUANTILE_H2G.out.quantile_results
        annotation_comparison = COMPARE_ANNOTATIONS.out.comparison_table
}