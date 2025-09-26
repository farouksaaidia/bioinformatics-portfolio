#!/bin/bash
# SnpEff variant annotation
# Input: raw VCF, genome database
# Output: annotated VCF

VCF=$1
DB=$2
OUT=$3

# Annotate variants
snpEff -v $DB $VCF > $OUT

# Optionally generate HTML summary
snpEff -v $DB $VCF -stats ${OUT%.vcf}_summary.html
