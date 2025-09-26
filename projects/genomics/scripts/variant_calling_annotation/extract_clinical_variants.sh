#!/bin/bash
# extract_clinical_variants.sh
# Input: annotated VCF
# Output: VCF containing only variants with clinical significance (ClinVar/COSMIC)

VCF=$1
OUT=$2

bcftools view -i 'INFO/CLNSIG!="."' $VCF -Ov -o $OUT
