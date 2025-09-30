#!/bin/bash
# subset_samples_vcf.sh → Extract specific samples from a multi-sample VCF
# Usage: ./subset_samples_vcf.sh input.vcf sample1,sample2 out.vcf
# Input: multi-sample VCF + comma-separated list of samples
# Output: subset VCF with only selected samples

VCF=$1
SAMPLES=$2
OUT=$3

bcftools view -s $SAMPLES -Ov -o $OUT $VCF
