#!/usr/bin/env bash
# Run samtools idxstats to get per-chromosome alignment stats

if [ $# -lt 1 ]; then
    echo "Usage: $0 <input.bam>"
    exit 1
fi

BAM=$1
OUT=${BAM%.bam}_idxstats.txt

samtools idxstats "$BAM" > "$OUT"

echo "Index stats report written to $OUT"
