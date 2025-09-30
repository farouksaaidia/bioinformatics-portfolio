#!/usr/bin/env bash
# Convert SAM to BAM
# Usage: ./sam_to_bam.sh input.sam output.bam

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input.sam output.bam"
    exit 1
fi

samtools view -S -b "$1" > "$2"
echo "✅ SAM converted to BAM: $2"
