#!/usr/bin/env bash
# Extract a genomic region from BAM or FASTA

if [ $# -lt 2 ]; then
    echo "Usage: $0 <input.bam/fa> <region: chr:start-end>"
    exit 1
fi

INPUT=$1
REGION=$2
OUT=${INPUT%.bam}_${REGION}.bam

if [[ "$INPUT" == *.bam ]]; then
    samtools view -b "$INPUT" "$REGION" > "$OUT"
    echo "Subset BAM written to $OUT"
elif [[ "$INPUT" == *.fa* ]]; then
    samtools faidx "$INPUT" "$REGION" > "${REGION}.fa"
    echo "Region FASTA written to ${REGION}.fa"
else
    echo "Unsupported input format. Must be BAM or FASTA."
    exit 1
fi
