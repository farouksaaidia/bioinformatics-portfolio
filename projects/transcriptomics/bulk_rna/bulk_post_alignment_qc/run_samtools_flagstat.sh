#!/bin/bash
# run_samtools_flagstat.sh
# Usage: ./run_samtools_flagstat.sh <bam_or_bam_dir> <outdir>
INPUT=$1
OUTDIR=${2:-flagstat_results}
mkdir -p "$OUTDIR"

if [ -z "$INPUT" ]; then
  echo "Usage: $0 <bam_or_bam_dir> <outdir>"
  exit 1
fi

if [ -d "$INPUT" ]; then
  for bam in "$INPUT"/*.bam; do
    [ -e "$bam" ] || continue
    sample=$(basename "$bam" .bam)
    samtools flagstat "$bam" > "$OUTDIR/${sample}_flagstat.txt"
  done
else
  sample=$(basename "$INPUT" .bam)
  samtools flagstat "$INPUT" > "$OUTDIR/${sample}_flagstat.txt"
fi
