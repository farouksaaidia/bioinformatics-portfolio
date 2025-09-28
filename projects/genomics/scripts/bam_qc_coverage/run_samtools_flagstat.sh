#!/usr/bin/env bash
# Run samtools flagstat to get basic alignment statistics

if [ $# -lt 1 ]; then
    echo "Usage: $0 <input.bam>"
    exit 1
fi

BAM=$1
OUT=${BAM%.bam}_flagstat.txt

samtools flagstat "$BAM" > "$OUT"

echo "Flagstat report written to $OUT"
