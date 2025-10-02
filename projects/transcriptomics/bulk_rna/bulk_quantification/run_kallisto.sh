#!/bin/bash
# Run Kallisto quantification for one or multiple FASTQ files
# Usage:
#   ./run_kallisto.sh sample_R1.fastq.gz sample_R2.fastq.gz index.idx output_dir
#   ./run_kallisto.sh fastq_folder index.idx output_dir

set -euo pipefail

INPUT=$1
INDEX=$2
OUTDIR=$3

mkdir -p "$OUTDIR"

if [ -d "$INPUT" ]; then
    echo "[INFO] Directory detected. Running Kallisto on all FASTQ pairs..."
    for R1 in "$INPUT"/*_R1.fastq.gz; do
        R2=${R1/_R1.fastq.gz/_R2.fastq.gz}
        SAMPLE=$(basename "$R1" _R1.fastq.gz)
        kallisto quant -i "$INDEX" -o "$OUTDIR/$SAMPLE" -b 100 -t 8 "$R1" "$R2"
    done
else
    R1=$INPUT
    R2=${R1/_R1.fastq.gz/_R2.fastq.gz}
    SAMPLE=$(basename "$R1" _R1.fastq.gz)
    kallisto quant -i "$INDEX" -o "$OUTDIR/$SAMPLE" -b 100 -t 8 "$R1" "$R2"
fi
