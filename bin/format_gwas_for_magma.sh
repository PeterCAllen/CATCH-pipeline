#!/usr/bin/env bash
set -euo pipefail

# Usage: format_gwas_for_magma.sh <GWAS_FILE> <BIM_FILE> <OUTPUT_PREFIX>
#   GWAS_FILE: path to input GWAS summary stats (gzipped)
#   BIM_FILE:  path to reference .bim file  
#   OUTPUT_PREFIX: prefix for output files
#
# This script handles GWAS files that either:
#   1. Already have rsIDs (e.g., rs123456) - uses them directly
#   2. Have positional IDs (e.g., 1:54490:G:A) - converts to rsIDs

GWAS_FILE=$1
BIM_FILE=$2
PREFIX=$3

echo "Formatting GWAS for MAGMA..." >&2

# First pass: detect if file has rsIDs by checking first SNP value
FIRST_SNP=$(gunzip -c ${GWAS_FILE} | awk '
NR==1 {
    # Find SNP column in header
    for(i=1; i<=NF; i++) {
        col_name = toupper($i);
        if(col_name == "SNP" || col_name == "RSID" || col_name == "RS" || col_name == "SNPID" || col_name == "MARKERNAME") {
            snp_col = i;
            break;
        }
    }
}
NR==2 {
    # Print first SNP value
    if(snp_col > 0) {
        print $snp_col;
    }
    exit;
}')

echo "  First SNP detected: ${FIRST_SNP}" >&2

# Check if conversion is needed
if [[ "$FIRST_SNP" =~ ^rs[0-9]+ ]]; then
    # Already has rsIDs - process directly
    echo "  ✓ GWAS already has rsIDs" >&2
    
    gunzip -c ${GWAS_FILE} | awk '
    NR==1 {
        # Detect column positions
        for(i=1; i<=NF; i++) {
            col_name = toupper($i);
            if(col_name == "CHR") chr_col = i;
            if(col_name == "BP" || col_name == "POS") bp_col = i;
            if(col_name == "SNP" || col_name == "RSID" || col_name == "RS" || col_name == "SNPID" || col_name == "MARKERNAME") snp_col = i;
            if(col_name == "P" || col_name == "PVAL" || col_name == "P_VALUE" || col_name == "PVALUE") p_col = i;
            if(col_name == "N" || col_name == "NEFF" || col_name == "N_TOTAL") n_col = i;
        }
        print "SNP", "P", "N";
        next;
    }
    NR>1 {
        if($p_col != "" && $p_col > 0 && $p_col <= 1 && $n_col != "") {
            print $snp_col, $p_col, $n_col;
        }
    }' > ${PREFIX}.pval
    
    gunzip -c ${GWAS_FILE} | awk '
    NR==1 {
        # Detect column positions
        for(i=1; i<=NF; i++) {
            col_name = toupper($i);
            if(col_name == "CHR") chr_col = i;
            if(col_name == "BP" || col_name == "POS") bp_col = i;
            if(col_name == "SNP" || col_name == "RSID" || col_name == "RS" || col_name == "SNPID" || col_name == "MARKERNAME") snp_col = i;
            if(col_name == "P" || col_name == "PVAL" || col_name == "P_VALUE" || col_name == "PVALUE") p_col = i;
        }
        next;
    }
    NR>1 {
        if($p_col != "" && $p_col > 0 && $p_col <= 1) {
            print $snp_col, $chr_col, $bp_col;
        }
    }' > ${PREFIX}.snp.loc
    
else
    # Needs conversion from positional IDs to rsIDs
    echo "  ✗ GWAS uses positional IDs - converting..." >&2
    
    # Create CHR:POS -> rsID mapping from reference BIM file
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
    NR==1 {
        # Detect column positions
        for(i=1; i<=NF; i++) {
            col_name = toupper($i);
            if(col_name == "CHR") chr_col = i;
            if(col_name == "BP" || col_name == "POS") bp_col = i;
            if(col_name == "P" || col_name == "PVAL" || col_name == "P_VALUE" || col_name == "PVALUE") p_col = i;
            if(col_name == "N" || col_name == "NEFF" || col_name == "N_TOTAL") n_col = i;
        }
        next;
    }
    NR>1 {
        pos_id = $chr_col":"$bp_col;
        if(pos_id in pos_rsid && $p_col != "" && $p_col > 0 && $p_col <= 1 && $n_col != "") {
            print pos_rsid[pos_id], $p_col, $n_col;
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
    NR==1 {
        # Detect column positions
        for(i=1; i<=NF; i++) {
            col_name = toupper($i);
            if(col_name == "CHR") chr_col = i;
            if(col_name == "BP" || col_name == "POS") bp_col = i;
            if(col_name == "P" || col_name == "PVAL" || col_name == "P_VALUE" || col_name == "PVALUE") p_col = i;
        }
        next;
    }
    NR>1 {
        pos_id = $chr_col":"$bp_col;
        if(pos_id in pos_rsid && $p_col != "" && $p_col > 0 && $p_col <= 1) {
            print pos_rsid[pos_id], $chr_col, $bp_col;
        }
    }' > ${PREFIX}.snp.loc
    
    rm -f pos_to_rsid.map
fi

echo "✓ Finished formatting GWAS for MAGMA. SNPs processed:" >&2
wc -l ${PREFIX}.pval >&2