#!/usr/bin/env bash
# Convert VCF to BCF or BCF to VCF
# Usage: ./convert_vcf_bcf.sh input.vcf output.bcf
#        ./convert_vcf_bcf.sh input.bcf output.vcf

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input.(vcf|bcf) output.(bcf|vcf)"
    exit 1
fi

bcftools view -O b -o "$2" "$1" 2>/dev/null || \
bcftools view -O v -o "$2" "$1"
echo "✅ Converted: $1 → $2"
