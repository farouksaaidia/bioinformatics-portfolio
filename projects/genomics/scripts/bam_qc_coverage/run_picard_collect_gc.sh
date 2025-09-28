#!/usr/bin/env bash
# Run Picard CollectGcBiasMetrics

if [ $# -lt 2 ]; then
    echo "Usage: $0 <input.bam> <reference.fa>"
    exit 1
fi

BAM=$1
REF=$2
OUT=${BAM%.bam}_gc_bias.txt
PDF=${BAM%.bam}_gc_bias.pdf
SUMMARY=${BAM%.bam}_gc_bias_summary.txt

picard CollectGcBiasMetrics \
    I="$BAM" \
    O="$OUT" \
    CHART="$PDF" \
    S="$SUMMARY" \
    R="$REF"

echo "GC bias metrics written to $OUT and $PDF"
