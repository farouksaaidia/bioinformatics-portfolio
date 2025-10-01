#!/bin/bash
# run_samtools_stats.sh
# Usage: ./run_samtools_stats.sh <bam_or_bam_dir> <outdir>
INPUT=$1
OUTDIR=${2:-stats_results}
mkdir -p "$OUTDIR"

if [ -z "$INPUT" ]; then
  echo "Usage: $0 <bam_or_bam_dir> <outdir>"
  exit 1
fi

if [ -d "$INPUT" ]; then
  for bam in "$INPUT"/*.bam; do
    [ -e "$bam" ] || continue
    sample=$(basename "$bam" .bam)
    samtools stats "$bam" > "$OUTDIR/${sample}_stats.txt"
  done
else
  sample=$(basename "$INPUT" .bam)
  samtools stats "$INPUT" > "$OUTDIR/${sample}_stats.txt"
fi
