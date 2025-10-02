#!/bin/bash
# Run Salmon quantification for one or multiple FASTQ files
# Usage:
#   ./run_salmon.sh sample_R1.fastq.gz sample_R2.fastq.gz index_dir output_dir
#   ./run_salmon.sh fastq_folder index_dir output_dir

set -euo pipefail

INPUT=$1
INDEX=$2
OUTDIR=$3

mkdir -p "$OUTDIR"

if [ -d "$INPUT" ]; then
    echo "[INFO] Directory detected. Running Salmon quant on all FASTQ pairs..."
    for R1 in "$INPUT"/*_R1.fastq.gz; do
        R2=${R1/_R1.fastq.gz/_R2.fastq.gz}
        SAMPLE=$(basename "$R1" _R1.fastq.gz)
        salmon quant -i "$INDEX" -l A -1 "$R1" -2 "$R2" -p 8 -o "$OUTDIR/$SAMPLE"
    done
else
    R1=$INPUT
    R2=${R1/_R1.fastq.gz/_R2.fastq.gz}
    SAMPLE=$(basename "$R1" _R1.fastq.gz)
    salmon quant -i "$INDEX" -l A -1 "$R1" -2 "$R2" -p 8 -o "$OUTDIR/$SAMPLE"
fi
