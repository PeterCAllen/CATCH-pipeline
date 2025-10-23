#!/usr/bin/env bash
set -euo pipefail

GWAS_FILE=$1
BIM_FILE=$2
PREFIX=$3

# Create a mapping file from reference panel: CHR:POS -> rsID
awk '{print $1":"$4, $2}' ${BIM_FILE} > pos_to_rsid.map

# Create p-value file with rsIDs
gunzip -c ${GWAS_FILE} | awk '
BEGIN {
    while((getline line < "pos_to_rsid.map") > 0) {
        split(line, arr, " ");
        pos_rsid[arr[1]] = arr[2];
    }
    close("pos_to_rsid.map");
    print "SNP", "P", "N";
}
NR>1 {
    pos_id = $1":"$2;
    if(pos_id in pos_rsid && $9 != "" && $9 > 0 && $9 <= 1) {
        print pos_rsid[pos_id], $9, $10;
    }
}' > ${PREFIX}.pval

# Create SNP location file with rsIDs
gunzip -c ${GWAS_FILE} | awk '
BEGIN {
    while((getline line < "pos_to_rsid.map") > 0) {
        split(line, arr, " ");
        pos_rsid[arr[1]] = arr[2];
    }
    close("pos_to_rsid.map");
}
NR>1 {
    pos_id = $1":"$2;
    if(pos_id in pos_rsid && $9 != "" && $9 > 0 && $9 <= 1) {
        print pos_rsid[pos_id], $1, $2;
    }
}' > ${PREFIX}.snp.loc

echo "Finished formatting GWAS for MAGMA. SNPs processed:"
wc -l ${PREFIX}.pval