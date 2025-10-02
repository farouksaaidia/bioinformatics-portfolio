#!/bin/bash
# Run HTSeq-count for one or multiple BAM files
# Usage:
#   ./run_htseq_count.sh input.bam annotation.gtf output_dir
#   ./run_htseq_count.sh /path/to/bam_folder annotation.gtf output_dir

set -euo pipefail

INPUT=$1
ANNOT=$2
OUTDIR=$3

mkdir -p "$OUTDIR"

if [ -d "$INPUT" ]; then
    echo "[INFO] Directory detected. Running on all BAM files in $INPUT"
    for BAM in "$INPUT"/*.bam; do
        SAMPLE=$(basename "$BAM" .bam)
        htseq-count -f bam -r pos -s no -a 10 "$BAM" "$ANNOT" > "$OUTDIR/${SAMPLE}_htseq_counts.txt"
    done
else
    SAMPLE=$(basename "$INPUT" .bam)
    htseq-count -f bam -r pos -s no -a 10 "$INPUT" "$ANNOT" > "$OUTDIR/${SAMPLE}_htseq_counts.txt"
fi
