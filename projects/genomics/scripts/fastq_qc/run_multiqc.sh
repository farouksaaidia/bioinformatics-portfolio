#!/bin/bash
# run_multiqc.sh → Aggregate FastQC results into a single report
# Usage: ./run_multiqc.sh [results_directory]

RESULTS_DIR=${1:-fastqc_results}

if [ ! -d "$RESULTS_DIR" ]; then
    echo "Error: directory $RESULTS_DIR not found!"
    exit 1
fi

echo "Running MultiQC on $RESULTS_DIR..."
multiqc "$RESULTS_DIR" -o multiqc_report

echo "MultiQC report generated in multiqc_report/"
