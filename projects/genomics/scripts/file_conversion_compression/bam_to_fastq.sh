#!/usr/bin/env bash
# Convert BAM to FASTQ
# Usage: ./bam_to_fastq.sh input.bam output_R1.fastq output_R2.fastq

set -euo pipefail

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 input.bam output_R1.fastq output_R2.fastq"
    exit 1
fi

samtools fastq -1 "$2" -2 "$3" "$1"
echo "✅ BAM converted to FASTQ: $2 and $3"
