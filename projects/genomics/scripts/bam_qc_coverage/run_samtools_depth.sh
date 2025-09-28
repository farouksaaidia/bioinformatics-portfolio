#!/usr/bin/env bash
# Run samtools depth to compute per-base coverage

if [ $# -lt 1 ]; then
    echo "Usage: $0 <input.bam>"
    exit 1
fi

BAM=$1
OUT=${BAM%.bam}_depth.txt

samtools depth -aa "$BAM" > "$OUT"

echo "Depth file written to $OUT"
