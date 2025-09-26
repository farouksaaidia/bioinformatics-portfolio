#!/bin/bash
# FreeBayes variant calling
# Input: sorted BAM, reference FASTA
# Output: raw VCF

BAM=$1
REF=$2
OUT=$3

# Call variants
freebayes -f $REF $BAM > $OUT
