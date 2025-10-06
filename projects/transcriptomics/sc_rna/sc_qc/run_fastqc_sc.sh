#!/usr/bin/env bash
# Purpose: Run FastQC on single or multiple scRNA-seq FASTQ samples

set -euo pipefail

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <output_dir> <fastq1> [fastq2 ...]"
    exit 1
fi

OUTDIR=$1
shift
mkdir -p "$OUTDIR"

for FASTQ in "$@"; do
    if [ ! -f "$FASTQ" ]; then
        echo "Error: File $FASTQ not found."
        exit 1
    fi
    fastqc -o "$OUTDIR" "$FASTQ"
done

echo "FastQC finished for all samples."
