#!/bin/bash
# Run featureCounts for one or multiple BAM files
# Usage:
#   ./run_featurecounts.sh input.bam annotation.gtf output_dir
#   ./run_featurecounts.sh /path/to/bam_folder annotation.gtf output_dir

set -euo pipefail

INPUT=$1
ANNOT=$2
OUTDIR=$3

mkdir -p "$OUTDIR"

if [ -d "$INPUT" ]; then
    echo "[INFO] Directory detected. Running on all BAM files in $INPUT"
    featureCounts -a "$ANNOT" -o "$OUTDIR/gene_counts.txt" "$INPUT"/*.bam
else
    echo "[INFO] Single file detected: $INPUT"
    SAMPLE=$(basename "$INPUT" .bam)
    featureCounts -a "$ANNOT" -o "$OUTDIR/${SAMPLE}_counts.txt" "$INPUT"
fi
