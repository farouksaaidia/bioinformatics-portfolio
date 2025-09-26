#!/bin/bash
# annotate_conservation_scores.sh
# Input: annotated VCF
# Output: VCF with conservation scores

VCF=$1
OUT=$2
DB="conservation_scores.vcf.gz"  # pre-indexed GERP/PhyloP VCF

bcftools annotate -a $DB -c INFO/GERP,INFO/PhyloP $VCF -Ov -o $OUT
