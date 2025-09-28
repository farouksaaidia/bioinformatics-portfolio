#!/usr/bin/env bash
# Run Picard CollectInsertSizeMetrics

if [ $# -lt 1 ]; then
    echo "Usage: $0 <input.bam>"
    exit 1
fi

BAM=$1
OUT=${BAM%.bam}_insert_size.txt
PDF=${BAM%.bam}_insert_size.pdf

picard CollectInsertSizeMetrics \
    I="$BAM" \
    O="$OUT" \
    H="$PDF" \
    M=0.5

echo "Insert size metrics written to $OUT and $PDF"
