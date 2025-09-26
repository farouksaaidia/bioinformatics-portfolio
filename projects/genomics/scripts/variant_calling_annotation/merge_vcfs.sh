#!/bin/bash
# Merge multiple VCF files
OUT=$1
shift
VCFS=$@

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <output.vcf> <vcf1> [vcf2 ...]"
    exit 1
fi

bcftools merge $VCFS -O v -o $OUT
echo "Merged VCF saved to $OUT"
