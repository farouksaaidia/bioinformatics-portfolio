#!/bin/bash
# FPKM normalization
# Input: counts.txt (gene_id, length, counts)
# Output: fpkm_matrix.txt

INPUT=$1
OUTPUT=$2

awk 'NR>1 { fpkm=($3*1e9)/($2*sum); print $1"\t"fpkm }' sum=$(awk 'NR>1 {sum+=$3} END {print sum}' $INPUT) $INPUT > $OUTPUT
