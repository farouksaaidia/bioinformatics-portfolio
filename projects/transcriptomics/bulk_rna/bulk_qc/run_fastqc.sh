#!/bin/bash
# Run FastQC on single or multiple FASTQ files

INPUT=$1
OUTPUT=${2:-fastqc_results}

mkdir -p "$OUTPUT"

if [ -d "$INPUT" ]; then
    # Process all FASTQ files in directory
    for fq in "$INPUT"/*.fastq*; do
        fastqc -o "$OUTPUT" "$fq"
    done
else
    # Process single FASTQ
    fastqc -o "$OUTPUT" "$INPUT"
fi
