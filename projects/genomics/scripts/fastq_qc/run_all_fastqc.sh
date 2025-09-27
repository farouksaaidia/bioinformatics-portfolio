#!/bin/bash
# run_all_fastqc.sh → Run FastQC on multiple FASTQ samples
# Usage: ./run_all_fastqc.sh /path/to/fastq_folder /path/to/output_folder

FASTQ_DIR=$1
OUT_DIR=$2

mkdir -p "$OUT_DIR"

for R1 in "$FASTQ_DIR"/*_1.fastq*; do
    R2="${R1/_1.fastq/_2.fastq}"
    SAMPLE=$(basename "$R1" "_1.fastq.gz")
    echo "Running FastQC for sample: $SAMPLE"
    fastqc -o "$OUT_DIR" -f fastq "$R1" "$R2"
done

echo "FastQC finished for all samples."
