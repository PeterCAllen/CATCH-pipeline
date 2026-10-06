#!/usr/bin/env bash
# Fetches all hg19/GRCh37 reference data the CATCH pipeline needs into data/reference/,
# matching the default paths in nextflow.config exactly. Run once, from a node with
# internet access (PBS compute nodes are typically air-gapped -- use a login node).
# Safe to re-run: skips any file/archive whose extracted contents are already present.
#
# Source: "CATCH Pipeline Reference" Zenodo record, DOI 10.5281/zenodo.23174917
# (CC-BY-4.0; see MANIFEST.md in that record for per-file provenance/citations -- none
# of this is original work, it's a repackaging of public reference data).
#
# Usage: bin/setup_reference_data.sh [target_dir]
#   target_dir defaults to ./data/reference (relative to the pipeline repo root).
set -euo pipefail

ZENODO_RECORD="23174917"
ZENODO_BASE="https://zenodo.org/api/records/${ZENODO_RECORD}/files"

TARGET_ROOT="${1:-data/reference}"
LDSC_DIR="${TARGET_ROOT}/LDSC_relevant_for_Peter"

mkdir -p "$LDSC_DIR"

# key: filename on Zenodo -> expected MD5 (from the record's own API metadata, checked
# 2026-10-06 -- re-verify with `curl -s https://zenodo.org/api/records/${ZENODO_RECORD}`
# if files are ever added/replaced on the record).
declare -A ZENODO_MD5=(
    [1000G_EUR_Phase3_plink.tar.gz]="dbeee1f5460046176f4a35d5770572bb"
    [1000G_Phase3_weights_hm3_no_MHC.tar]="11dffb2e860906c58f1af8cd3941217e"
    [baseline_v1.1_thin_annot.tar]="000ee1abbb49f627e5bee0a43fc66489"
    [cellect_ldsc_support_files.tar.gz]="6546f2e91eddc0446c3a9d7df4c08b3a"
    [GRCh37pos.unique.ensgid.56778.cellex.format.txt]="887833ed196e601fd117fec279f9b904"
    [geneMatrix.tsv.gz]="66cd08848eef5256071875628a4dbced"
    [g1000_eur_magma.tar.gz]="7219d0099e9bd07ce20b1ffeaa96c701"
)

# Downloads $1 from Zenodo to $2 and checksum-verifies it, unless $2 already exists.
fetch_and_verify() {
    local fname="$1" dest="$2"
    if [ -s "$dest" ]; then
        echo "  already present: $dest (skipping download)"
        return
    fi
    echo "  downloading $fname -> $dest"
    curl -sL -o "$dest" "${ZENODO_BASE}/${fname}/content"

    local expected="${ZENODO_MD5[$fname]:-}"
    if [ -n "$expected" ]; then
        local actual
        actual=$(md5sum "$dest" | cut -d' ' -f1)
        if [ "$actual" != "$expected" ]; then
            echo "ERROR: checksum mismatch for $fname (expected $expected, got $actual)" >&2
            rm -f "$dest"
            exit 1
        fi
    fi
}

# Fetches archive $1 and extracts it into $2 with tar flags $3 (e.g. "-xzf" or "-xf"),
# unless sentinel file $4 (a file expected to exist post-extraction) is already present.
fetch_extract_and_verify() {
    local fname="$1" extract_dir="$2" tar_flags="$3" sentinel="$4"
    if [ -s "$sentinel" ]; then
        echo "  already extracted: $sentinel (skipping)"
        return
    fi
    local tmp="${extract_dir}/${fname}"
    fetch_and_verify "$fname" "$tmp"
    mkdir -p "$extract_dir"
    tar "$tar_flags" "$tmp" -C "$extract_dir"
    rm "$tmp"
}

echo "=== 1000G EUR Phase 3 PLINK panel (LDSC/conLDSC) ==="
fetch_extract_and_verify "1000G_EUR_Phase3_plink.tar.gz" "$LDSC_DIR" "-xzf" \
    "${LDSC_DIR}/1000G_EUR_Phase3_plink/1000G.EUR.QC.1.bed"

echo "=== 1000G Phase 3 HM3 LD-score weights ==="
fetch_extract_and_verify "1000G_Phase3_weights_hm3_no_MHC.tar" "$LDSC_DIR" "-xf" \
    "${LDSC_DIR}/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.1.l2.ldscore.gz"

echo "=== baseline v1.1 thin-annot LD score model ==="
fetch_extract_and_verify "baseline_v1.1_thin_annot.tar" "$LDSC_DIR" "-xf" \
    "${LDSC_DIR}/baseline_v1.1_thin_annot/baseline.1.l2.ldscore.gz"

echo "=== CELLECT-LDSC support files (print_snps.txt, GRCh37-chr-sizes.txt, w_hm3.snplist) ==="
fetch_extract_and_verify "cellect_ldsc_support_files.tar.gz" "$LDSC_DIR" "-xzf" \
    "${LDSC_DIR}/print_snps.txt"

echo "=== conLDSC gene coordinate reference (GRCh37pos...56778...) ==="
fetch_and_verify "GRCh37pos.unique.ensgid.56778.cellex.format.txt" \
    "${LDSC_DIR}/GRCh37pos.unique.ensgid.56778.cellex.format.txt"

echo "=== geneMatrix.tsv.gz (required for every run -- feeds PREPARE_GENE_COORDS) ==="
fetch_and_verify "geneMatrix.tsv.gz" "${TARGET_ROOT}/geneMatrix.tsv.gz"

echo "=== MAGMA g1000_eur panel (mBAT-combo LD reference) ==="
fetch_extract_and_verify "g1000_eur_magma.tar.gz" "$TARGET_ROOT" "-xzf" \
    "${TARGET_ROOT}/g1000_eur_magma/g1000_eur.bed"

echo
echo "=== Verifying everything nextflow.config's defaults expect is present ==="
missing=0
check() {
    if [ ! -s "$1" ]; then
        echo "  MISSING: $1" >&2
        missing=1
    else
        echo "  OK: $1"
    fi
}
check "${TARGET_ROOT}/geneMatrix.tsv.gz"
check "${LDSC_DIR}/1000G_EUR_Phase3_plink/1000G.EUR.QC.1.bed"
check "${LDSC_DIR}/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.1.l2.ldscore.gz"
check "${LDSC_DIR}/baseline_v1.1_thin_annot/baseline.1.l2.ldscore.gz"
check "${LDSC_DIR}/print_snps.txt"
check "${LDSC_DIR}/GRCh37-chr-sizes.txt"
check "${LDSC_DIR}/w_hm3.snplist"
check "${LDSC_DIR}/GRCh37pos.unique.ensgid.56778.cellex.format.txt"
check "${TARGET_ROOT}/g1000_eur_magma/g1000_eur.bed"

if [ "$missing" -ne 0 ]; then
    echo
    echo "One or more reference files are missing -- see above. The pipeline will fail" >&2
    echo "at the corresponding step until these are resolved." >&2
    exit 1
fi

echo
echo "All reference data present under ${TARGET_ROOT}/. Next: build the containers"
echo "(see README.md), then run the pipeline."
