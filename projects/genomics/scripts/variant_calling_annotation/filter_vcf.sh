#!/bin/bash
# Filter VCF based on quality, depth, and allele frequency
VCF=$1
OUT=$2

if [[ ! -f $VCF ]]; then
    echo "Error: $VCF not found."
    exit 1
fi

# Example filters: QUAL>20, DP>10, AF>0.05
bcftools filter -i 'QUAL>20 && DP>10 && AF>0.05' $VCF -o $OUT
echo "Filtered VCF saved to $OUT"
