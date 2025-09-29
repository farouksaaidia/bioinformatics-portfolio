#!/usr/bin/env bash
# Extract random reads from FASTQ

if [ $# -lt 2 ]; then
    echo "Usage: $0 <input.fastq.gz> <num_reads>"
    exit 1
fi

FASTQ=$1
NUM=$2
OUT=${FASTQ%.fastq.gz}_subset_${NUM}.fastq.gz

zcat "$FASTQ" | paste - - - - | shuf -n "$NUM" | tr '\t' '\n' | gzip > "$OUT"

echo "Random subset FASTQ written to $OUT"
