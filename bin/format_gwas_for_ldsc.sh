#!/usr/bin/env bash
# set -euo pipefail

# Usage: format_gwas_for_ldsc.sh <GWAS_FILE> <BIM_FILE> <OUTPUT_PREFIX>
#   GWAS_FILE: path to input GWAS summary stats (gzipped)
#   BIM_FILE:  path to reference .bim file
#   OUTPUT_PREFIX: prefix for output files

GWAS_FILE=$1
BIM_FILE=$2
PREFIX=$3

echo "========================================" >&2
echo "Formatting GWAS for LDSC" >&2
echo "========================================" >&2
echo "Input: ${GWAS_FILE}" >&2
echo "Reference: ${BIM_FILE}" >&2
echo "" >&2

# First pass: detect if file has rsIDs by checking first SNP value
echo "Analyzing GWAS file format..." >&2
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

# Check if SNP column already contains rsIDs
if [[ "$FIRST_SNP" =~ ^rs[0-9]+ ]]; then
    # SNP column already has rsIDs - just copy the file as-is
    echo "✓ GWAS file already contains rsIDs" >&2
    echo "  No conversion needed - passing through..." >&2
    
    # Convert tabs to spaces for consistency with converted output
    gunzip -c ${GWAS_FILE} | awk '{$1=$1; print}' > ${PREFIX}_with_rsids.txt
    
else
    # SNP column has positional IDs - convert to rsIDs
    echo "✗ GWAS file uses positional IDs" >&2
    echo "  Converting to rsIDs from reference panel..." >&2
    echo "" >&2
    
    # Create a mapping file from reference panel: CHR:POS -> rsID
    awk '{print $1":"$4, $2}' ${BIM_FILE} > pos_to_rsid.map
    
    # Convert the GWAS file to use rsIDs instead of positional IDs
    # This awk script auto-detects columns and converts in one pass
    gunzip -c ${GWAS_FILE} | awk '
    BEGIN {
        # Load position to rsID mapping
        while((getline line < "pos_to_rsid.map") > 0) {
            split(line, arr, " ");
            pos_rsid[arr[1]] = arr[2];
        }
        close("pos_to_rsid.map");
        converted = 0;
        not_found = 0;
    }
    NR==1 {
        # Detect column positions from header
        for(i=1; i<=NF; i++) {
            col_name = toupper($i);
            if(col_name == "CHR") chr_col = i;
            if(col_name == "BP" || col_name == "POS") bp_col = i;
            if(col_name == "SNP" || col_name == "RSID" || col_name == "RS" || col_name == "SNPID" || col_name == "MARKERNAME") snp_col = i;
        }
        # Print header as-is
        print;
        next;
    }
    NR>1 {
        # Create position key from CHR and BP columns
        pos_id = $chr_col":"$bp_col;
        
        # Check if position exists in reference panel
        if(pos_id in pos_rsid) {
            # Replace SNP column with rsID
            $snp_col = pos_rsid[pos_id];
            converted++;
        } else {
            # Keep original if no match found
            not_found++;
        }
        print;
    }
    END {
        print "  ✓ Converted " converted " variants to rsIDs" > "/dev/stderr";
        if(not_found > 0) {
            print "  ⚠ " not_found " variants not found in reference" > "/dev/stderr";
        }
    }' > ${PREFIX}_with_rsids.txt
    
    rm -f pos_to_rsid.map
fi

echo "" >&2
echo "✓ Finished formatting GWAS for LDSC" >&2
echo "  Output: ${PREFIX}_with_rsids.txt" >&2
echo "" >&2

# Show a preview of the output
echo "Preview (first 3 rows, first 6 columns):" >&2
head -n 3 ${PREFIX}_with_rsids.txt | cut -f1-6 >&2
