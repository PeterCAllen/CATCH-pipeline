#!/usr/bin/env bash
set -euo pipefail

# Usage: format_gwas_for_ldsc.sh <GWAS_FILE> <BIM_FILE> <OUTPUT_PREFIX>
#   GWAS_FILE: path to input GWAS summary stats (gzipped)
#   BIM_FILE:  path to reference .bim file
#   OUTPUT_PREFIX: prefix for output files

GWAS_FILE=$1
BIM_FILE=$2
PREFIX=$3

# Create a mapping file from reference panel: CHR:POS -> rsID
awk '{print $1":"$4, $2}' ${BIM_FILE} > pos_to_rsid.map

# Convert the GWAS file to use rsIDs instead of positional IDs
gunzip -c ${GWAS_FILE} | awk '
BEGIN {
	# Load position to rsID mapping
	while((getline line < "pos_to_rsid.map") > 0) {
		split(line, arr, " ");
		pos_rsid[arr[1]] = arr[2];
	}
	close("pos_to_rsid.map");
}
NR==1 {
	# Print header with SNP as first column, rename original positional SNP column
	print "SNP", $2, "SNP_POS", $4, $5, $6, $7, $8, $9, $10, $11;
}
NR>1 {
	pos_id = $1":"$2;
	# Check if position exists in reference panel
	if(pos_id in pos_rsid) {
		# Print rsID, BP, original positional ID (renamed), and remaining columns
		print pos_rsid[pos_id], $2, $1":"$2, $4, $5, $6, $7, $8, $9, $10, $11;
	}
}' > ${PREFIX}_with_rsids.txt

echo "Finished formatting GWAS for LDSC. Output: ${PREFIX}_with_rsids.txt"
