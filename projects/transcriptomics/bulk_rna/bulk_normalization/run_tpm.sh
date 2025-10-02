#!/bin/bash
# Simple TPM normalization from counts.txt
# Input: counts.txt (tab-delimited: gene_id, length, counts)
# Output: tpm_matrix.txt

INPUT=$1
OUTPUT=$2

awk 'NR>1 { rpm=$3*1e6/sum; tpm=rpm/($2/1000); print $1"\t"tpm }' sum=$(awk 'NR>1 {sum+=$3} END {print sum}' $INPUT) $INPUT > $OUTPUT
