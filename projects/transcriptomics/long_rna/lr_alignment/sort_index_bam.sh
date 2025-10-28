#!/usr/bin/env bash
set -euo pipefail

INPUT=$1
THREADS=${2:-8}

[ ! -f "$INPUT" ] && { echo "Error: BAM file missing"; exit 1; }

samtools sort -@ "$THREADS" "$INPUT" -o "${INPUT%.bam}.sorted.bam"
samtools index "${INPUT%.bam}.sorted.bam"

echo "✅ Sorted and indexed → ${INPUT%.bam}.sorted.bam"
