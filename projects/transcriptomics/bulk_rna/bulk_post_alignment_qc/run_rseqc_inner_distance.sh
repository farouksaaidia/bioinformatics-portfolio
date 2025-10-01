#!/bin/bash
# run_rseqc_inner_distance.sh
# Usage: ./run_rseqc_inner_distance.sh <bam_or_bam_dir> <gtf> <outdir>
INPUT=$1
GTF=$2
OUTDIR=${3:-inner_distance_results}
mkdir -p "$OUTDIR"

if [ -z "$INPUT" ] || [ -z "$GTF" ]; then
  echo "Usage: $0 <bam_or_bam_dir> <gtf> <outdir>"
  exit 1
fi

if [ -d "$INPUT" ]; then
  for bam in "$INPUT"/*.bam; do
    [ -e "$bam" ] || continue
    sample=$(basename "$bam" .bam)
    inner_distance.py -r "$GTF" -i "$bam" -o "$OUTDIR/${sample}"
  done
else
  sample=$(basename "$INPUT" .bam)
  inner_distance.py -r "$GTF" -i "$INPUT" -o "$OUTDIR/${sample}"
fi
