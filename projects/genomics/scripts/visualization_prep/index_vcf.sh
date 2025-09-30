#!/bin/bash
# index_vcf.sh → Compress and index VCF for visualization
# Input: VCF
# Output: bgzipped VCF + .tbi index

VCF=$1

bgzip -c $VCF > ${VCF}.gz
tabix -p vcf ${VCF}.gz
