#!/bin/bash
# run_rseqc_infer_experiment.sh
# Usage: ./run_rseqc_infer_experiment.sh <bam_or_bam_dir> <gtf> <outdir>
INPUT=$1
GTF=$2
OUTDIR=${3:-infer_exp_results}
mkdir -p "$OUTDIR"

if [ -z "$INPUT" ] || [ -z "$GTF" ]; then
  echo "Usage: $0 <bam_or_bam_dir> <gtf> <outdir>"
  exit 1
fi

if [ -d "$INPUT" ]; then
  for bam in "$INPUT"/*.bam; do
    [ -e "$bam" ] || continue
    sample=$(basename "$bam" .bam)
    infer_experiment.py -r "$GTF" -i "$bam" > "$OUTDIR/${sample}_infer.txt"
  done
else
  sample=$(basename "$INPUT" .bam)
  infer_experiment.py -r "$GTF" -i "$INPUT" > "$OUTDIR/${sample}_infer.txt"
fi
