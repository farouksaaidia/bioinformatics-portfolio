#!/bin/bash
# run_rseqc_geneBody.sh
# Usage: ./run_rseqc_geneBody.sh <bam_or_bam_dir> <gtf> <outdir>
INPUT=$1
GTF=$2
OUTDIR=${3:-genebody_results}
mkdir -p "$OUTDIR"

if [ -z "$INPUT" ] || [ -z "$GTF" ]; then
  echo "Usage: $0 <bam_or_bam_dir> <gtf> <outdir>"
  exit 1
fi

if [ -d "$INPUT" ]; then
  for bam in "$INPUT"/*.bam; do
    [ -e "$bam" ] || continue
    sample=$(basename "$bam" .bam)
    geneBody_coverage.py -r "$GTF" -i "$bam" -o "$OUTDIR/${sample}"
  done
else
  sample=$(basename "$INPUT" .bam)
  geneBody_coverage.py -r "$GTF" -i "$INPUT" -o "$OUTDIR/${sample}"
fi
