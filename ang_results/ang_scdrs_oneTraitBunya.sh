#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --job-name=scDRS_one_trait_bunya
#SBATCH --partition=general
#SBATCH --account=a_imb_cpdg
#SBATCH --error=slurm_one_trait_%j.error
#SBATCH --output=slurm_one_trait_%j.out

set -euo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
CONDA_ENV_PATH="${CONDA_ENV_PATH:-/scratch/user/uqali4/micromamba/envs/scRNA_env}"
conda activate "${CONDA_ENV_PATH}"

PIPELINE="${PIPELINE:-A}"
DATASET_NAME="${DATASET_NAME:-}"
TRAIT="${TRAIT:-BMI}"
METHOD="${METHOD:-MAGMA}"

BASE_SCDATA="${BASE_SCDATA:-/QRISdata/Q9052/scData/preprocess}"
BASE_RESULTS="${BASE_RESULTS:-/QRISdata/Q9052/scDRS_res}"
BASE_GS_ROOT="${BASE_GS_ROOT:-/QRISdata/Q9052/benchmark}"
FILE_SUFFIX="${FILE_SUFFIX:-_10kb}"

if [[ -z "${DATASET_NAME}" ]]; then
  echo "ERROR: DATASET_NAME is required"
  exit 1
fi

case "${METHOD}" in
  MAGMA) GS_DIR="${BASE_GS_ROOT}/MAGMA_1kg_eur_benchmark" ; SUFFIX="" ;;
  mBAT)  GS_DIR="${BASE_GS_ROOT}/mBAT_1kg_eur_benchmark" ; SUFFIX="_mBAT" ;;
  *) echo "ERROR: METHOD must be MAGMA or mBAT"; exit 1 ;;
esac

case "${PIPELINE}" in
  A)
    score_folder="${BASE_RESULTS}/A_baseline_step1${SUFFIX}/${DATASET_NAME}"
    downstream_folder="${BASE_RESULTS}/A_baseline_step2${SUFFIX}/${DATASET_NAME}"
    ;;
  B)
    score_folder="${BASE_RESULTS}/B_expr_universe_step1${SUFFIX}/${DATASET_NAME}"
    downstream_folder="${BASE_RESULTS}/B_expr_universe_step2${SUFFIX}/${DATASET_NAME}"
    ;;
  C)
    score_folder="${BASE_RESULTS}/C_gwas_universe_step1${SUFFIX}/${DATASET_NAME}"
    downstream_folder="${BASE_RESULTS}/C_gwas_universe_step2${SUFFIX}/${DATASET_NAME}"
    ;;
  *)
    echo "ERROR: PIPELINE must be A, B, or C"
    exit 1
    ;;
esac

h5ad_file="${BASE_SCDATA}/${DATASET_NAME}.h5ad"
cov_file="${BASE_SCDATA}/${DATASET_NAME}.cov"
gs_file="${GS_DIR}/${DATASET_NAME}_${TRAIT}${FILE_SUFFIX}.gs"

case "${DATASET_NAME}" in
  2019_Smillie_normal_min20|2019_Smillie_normal_min20_pc)
    GROUP_ANALYSIS="cell_type"; FLAG_RAW_COUNT="False" ;;
  CARE_snRNA_Heart|CARE_snRNA_Heart_pc)
    GROUP_ANALYSIS="celltype"; FLAG_RAW_COUNT="False" ;;
  Cheng_2018_Cell_Reports_updated|Cheng_2018_Cell_Reports_pc)
    GROUP_ANALYSIS="cell_type"; FLAG_RAW_COUNT="False" ;;
  Fasolino_2022_Nat_Metab_normal_only|Fasolino_2022_Nat_Metab_normal_only_pc)
    GROUP_ANALYSIS="cell_type"; FLAG_RAW_COUNT="False" ;;
  FBM_Jardine_2021_Nature|FBM_Jardine_2021_Nature_pc)
    GROUP_ANALYSIS="cell_type"; FLAG_RAW_COUNT="False" ;;
  HLCA_core_healthy_LP|HLCA_core_healthy_LP_pc)
    GROUP_ANALYSIS="cell_type"; FLAG_RAW_COUNT="True" ;;
  human_liver_atlas_Guilliams_2022_cell|human_liver_atlas_Guilliams_2022_cell_pc)
    GROUP_ANALYSIS="cell_type"; FLAG_RAW_COUNT="False" ;;
  Kamath_2022_normal|Kamath_2022_normal_pc)
    GROUP_ANALYSIS="Cell_type_refined"; FLAG_RAW_COUNT="True" ;;
  TabulaSapiens_pc_minCell20|TabulaSapiens_pc_1to1_TMS_min20)
    GROUP_ANALYSIS="cell_ontology_class"; FLAG_RAW_COUNT="True" ;;
  TMS_pc_1to1_TS_min20)
    GROUP_ANALYSIS="cell_ontology_class"; FLAG_RAW_COUNT="True" ;;
  *)
    echo "ERROR: Unknown DATASET_NAME=${DATASET_NAME}"
    exit 1
    ;;
esac

mkdir -p "${score_folder}" "${downstream_folder}"

echo "============================================"
echo "Bunya | METHOD=${METHOD} PIPELINE=${PIPELINE} DATASET=${DATASET_NAME} TRAIT=${TRAIT}"
echo "============================================"

[[ -f "${gs_file}" ]] || { echo "ERROR: missing GS file: ${gs_file}"; exit 1; }
[[ -f "${h5ad_file}" ]] || { echo "ERROR: missing h5ad: ${h5ad_file}"; exit 1; }
[[ -f "${cov_file}" ]] || { echo "ERROR: missing cov: ${cov_file}"; exit 1; }

OUT_TEMP="${score_folder}/.tmp_${TRAIT}_${METHOD}"
rm -rf "${OUT_TEMP}"
mkdir -p "${OUT_TEMP}"

scdrs compute-score \
  --h5ad-file "${h5ad_file}" \
  --h5ad-species human \
  --gs-file "${gs_file}" \
  --gs-species human \
  --out-folder "${OUT_TEMP}" \
  --cov-file "${cov_file}" \
  --flag-filter-data True \
  --flag-raw-count "${FLAG_RAW_COUNT}" \
  --n-ctrl 1000 \
  --flag-return-ctrl-raw-score False \
  --flag-return-ctrl-norm-score True

final_full="${score_folder}/${TRAIT}${FILE_SUFFIX}.full_score.gz"
final_score="${score_folder}/${TRAIT}${FILE_SUFFIX}.score.gz"
scored=$(find "${OUT_TEMP}" -maxdepth 1 -name "*.full_score.gz" | head -n 1)
[[ -n "${scored}" ]] && mv "${scored}" "${final_full}"
scored2=$(find "${OUT_TEMP}" -maxdepth 1 -name "*.score.gz" | head -n 1)
[[ -n "${scored2}" ]] && mv "${scored2}" "${final_score}"
rm -rf "${OUT_TEMP}"

[[ -f "${final_full}" ]] || { echo "ERROR: no final full_score produced"; exit 1; }

scdrs perform-downstream \
  --h5ad-file "${h5ad_file}" \
  --score-file "${final_full}" \
  --out-folder "${downstream_folder}" \
  --group-analysis "${GROUP_ANALYSIS}" \
  --gene-analysis \
  --flag-filter-data True \
  --flag-raw-count "${FLAG_RAW_COUNT}"

echo "Done: ${METHOD} | ${PIPELINE} | ${DATASET_NAME} | ${TRAIT}"
