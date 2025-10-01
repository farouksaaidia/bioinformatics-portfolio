#!/bin/bash
# run_bamqc.sh
# Usage: ./run_bamqc.sh <bam_or_bam_dir> <outdir>
INPUT=$1
OUTDIR=${2:-bamqc_results}
mkdir -p "$OUTDIR"

if [ -z "$INPUT" ]; then
  echo "Usage: $0 <bam_or_bam_dir> <outdir>"
  exit 1
fi

if [ -d "$INPUT" ]; then
  for bam in "$INPUT"/*.bam; do
    [ -e "$bam" ] || continue
    bamqc -bam "$bam" -outdir "$OUTDIR"
  done
else
  bamqc -bam "$INPUT" -outdir "$OUTDIR"
fi
