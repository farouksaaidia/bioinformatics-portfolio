#!/bin/bash
# run_picard_markdups.sh
# Usage: ./run_picard_markdups.sh <bam_or_bam_dir> <outdir>
INPUT=$1
OUTDIR=${2:-markdups_results}
mkdir -p "$OUTDIR"

if [ -z "$INPUT" ]; then
  echo "Usage: $0 <bam_or_bam_dir> <outdir>"
  exit 1
fi

if [ -d "$INPUT" ]; then
  for bam in "$INPUT"/*.bam; do
    [ -e "$bam" ] || continue
    sample=$(basename "$bam" .bam)
    picard MarkDuplicates I="$bam" O="$OUTDIR/${sample}_dedup.bam" \
      M="$OUTDIR/${sample}_dup_metrics.txt" REMOVE_DUPLICATES=false CREATE_INDEX=true
  done
else
  sample=$(basename "$INPUT" .bam)
  picard MarkDuplicates I="$INPUT" O="$OUTDIR/${sample}_dedup.bam" \
    M="$OUTDIR/${sample}_dup_metrics.txt" REMOVE_DUPLICATES=false CREATE_INDEX=true
fi
