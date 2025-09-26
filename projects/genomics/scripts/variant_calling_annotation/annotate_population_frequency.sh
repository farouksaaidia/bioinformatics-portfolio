#!/bin/bash
# annotate_population_frequency.sh
# Input: annotated VCF
# Output: VCF with population frequency annotations

VCF=$1
OUT=$2
DB="gnomad.vcf.gz"  # path to pre-indexed gnomAD VCF

bcftools annotate -a $DB -c INFO/AF $VCF -Ov -o $OUT
