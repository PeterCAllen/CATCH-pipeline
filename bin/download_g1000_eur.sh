#!/usr/bin/env bash
# Downloads MAGMA's g1000_eur reference panel (503 EUR samples, ~22.6M SNPs, one
# combined whole-genome PLINK bed/bim/fam) for the mBAT-combo branch.
#
# Run this once from a node with internet access (e.g. a login node -- PBS compute
# nodes are typically air-gapped) before running the pipeline with --run_scdrs or
# --run_seismic. Output lands at data/reference/g1000_eur_magma/g1000_eur.{bed,bim,fam},
# matching the default --ref_hg19_mbat_plink_prefix in nextflow.config.
#
# Source: https://cncr.nl/research/magma/ ("European reference data (1000 Genomes Phase 3)")
set -euo pipefail

OUT_DIR="${1:-data/reference/g1000_eur_magma}"
URL="https://vu.data.surf.nl/index.php/s/VZNByNwpD8qqINe/download"

mkdir -p "$OUT_DIR"
tmp_zip="$(mktemp --suffix=.zip)"
trap 'rm -f "$tmp_zip"' EXIT

echo "Downloading g1000_eur.zip (~488MB) to $tmp_zip ..."
curl -sL -o "$tmp_zip" "$URL"

echo "Extracting into $OUT_DIR ..."
unzip -o "$tmp_zip" -d "$OUT_DIR"

for ext in bed bim fam; do
    if [ ! -s "$OUT_DIR/g1000_eur.$ext" ]; then
        echo "ERROR: $OUT_DIR/g1000_eur.$ext missing or empty after extraction" >&2
        exit 1
    fi
done

echo "Done. $(wc -l < "$OUT_DIR/g1000_eur.fam") samples, $(wc -l < "$OUT_DIR/g1000_eur.bim") SNPs."
