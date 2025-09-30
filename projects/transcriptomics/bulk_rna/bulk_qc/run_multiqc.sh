#!/bin/bash
# Aggregate FastQC or other QC reports with MultiQC

INPUT_DIR=$1
OUTPUT_DIR=${2:-multiqc_report}

mkdir -p "$OUTPUT_DIR"
multiqc "$INPUT_DIR" -o "$OUTPUT_DIR"
