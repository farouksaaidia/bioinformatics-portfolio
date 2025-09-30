#!/bin/bash
# subset_variants_vcf.sh → Keep/exclude specific variants or regions from a VCF
# Usage: ./subset_variants_vcf.sh input.vcf regions.bed out.vcf
# Input: VCF + BED/region file
# Output: VCF containing only variants in the specified regions

VCF=$1
REGIONS=$2
OUT=$3

bcftools view -R $REGIONS -Ov -o $OUT $VCF
