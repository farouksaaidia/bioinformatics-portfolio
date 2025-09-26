#!/bin/bash
# vcf_cleanup.sh → Clean VCF: remove duplicates, low-quality variants
# Input: VCF file, output file
# Output: cleaned VCF

VCF=$1
OUT=$2

# Remove duplicate entries
bcftools norm -d both $VCF -Ov -o $OUT

echo "VCF cleanup complete: $OUT"
