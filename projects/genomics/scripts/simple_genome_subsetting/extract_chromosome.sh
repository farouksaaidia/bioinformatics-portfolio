#!/usr/bin/env bash
# Extract a single chromosome from BAM or FASTA

if [ $# -lt 2 ]; then
    echo "Usage: $0 <input.bam/fa> <chromosome>"
    exit 1
fi

INPUT=$1
CHROM=$2
OUT=${INPUT%.bam}_${CHROM}.bam

if [[ "$INPUT" == *.bam ]]; then
    samtools view -b "$INPUT" "$CHROM" > "$OUT"
    echo "Subset BAM written to $OUT"
elif [[ "$INPUT" == *.fa* ]]; then
    samtools faidx "$INPUT" "$CHROM" > "${CHROM}.fa"
    echo "Chromosome FASTA written to ${CHROM}.fa"
else
    echo "Unsupported input format. Must be BAM or FASTA."
    exit 1
fi
