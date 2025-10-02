#!/bin/bash
# Counts Per Million normalization
# Input: counts.txt (gene_id, counts)
# Output: cpm_matrix.txt

INPUT=$1
OUTPUT=$2

awk 'NR>1 { cpm=($2*1e6)/sum; print $1"\t"cpm }' sum=$(awk 'NR>1 {sum+=$2} END {print sum}' $INPUT) $INPUT > $OUTPUT
