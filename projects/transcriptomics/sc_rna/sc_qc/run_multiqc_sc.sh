#!/usr/bin/env bash
# Purpose: Aggregate FastQC reports into a single MultiQC report

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_dir_of_fastqc_results> <output_dir>"
    exit 1
fi

INPUT_DIR=$1
OUTPUT_DIR=$2
mkdir -p "$OUTPUT_DIR"

multiqc "$INPUT_DIR" -o "$OUTPUT_DIR"

echo "MultiQC report generated at $OUTPUT_DIR"
