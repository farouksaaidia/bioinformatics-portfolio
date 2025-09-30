#!/usr/bin/env bash
# Convert BAM to SAM
# Usage: ./bam_to_sam.sh input.bam output.sam

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input.bam output.sam"
    exit 1
fi

samtools view -h "$1" > "$2"
echo "✅ BAM converted to SAM: $2"
