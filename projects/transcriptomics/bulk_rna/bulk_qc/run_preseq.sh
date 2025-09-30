#!/bin/bash
# Run Preseq to estimate library complexity from BAM file(s)

INPUT=$1
OUTPUT=${2:-preseq_results}

mkdir -p "$OUTPUT"

if [ -d "$INPUT" ]; then
    for bam in "$INPUT"/*.bam; do
        prefix=$(basename "$bam" .bam)
        preseq lc_extrap -B "$bam" -o "$OUTPUT/${prefix}_preseq.txt"
    done
else
    prefix=$(basename "$INPUT" .bam)
    preseq lc_extrap -B "$INPUT" -o "$OUTPUT/${prefix}_preseq.txt"
fi
