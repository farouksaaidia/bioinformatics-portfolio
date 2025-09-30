#!/usr/bin/env bash
# Sort and index a BAM file
# Usage: ./bam_sort_index.sh input.bam output.sorted.bam

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input.bam output.sorted.bam"
    exit 1
fi

samtools sort -o "$2" "$1"
samtools index "$2"

echo "✅ BAM sorted and indexed: $2 and $2.bai"
