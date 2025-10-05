#!/bin/bash
# =========================================
# run_alevin.sh
# Preprocessing pipeline for scRNA-seq using Salmon Alevin
# =========================================

set -euo pipefail

usage() {
    echo "Usage: $0 -i <sample_id> -1 <R1.fastq.gz> -2 <R2.fastq.gz> -r <salmon_index> -o <output_dir>"
    exit 1
}

while getopts ":i:1:2:r:o:h" opt; do
    case $opt in
        i) SAMPLE=$OPTARG ;;
        1) R1=$OPTARG ;;
        2) R2=$OPTARG ;;
        r) INDEX=$OPTARG ;;
        o) OUTDIR=$OPTARG ;;
        h) usage ;;
    esac
done

if [[ -z "${SAMPLE:-}" || -z "${R1:-}" || -z "${R2:-}" || -z "${INDEX:-}" || -z "${OUTDIR:-}" ]]; then
    echo "[ERROR] Missing arguments"; usage
fi

echo "[$(date)] Running Salmon Alevin for $SAMPLE"
salmon alevin -l ISR -1 "$R1" -2 "$R2" --chromium -i "$INDEX" -p 8 -o "$OUTDIR/$SAMPLE"
echo "[$(date)] Finished Salmon Alevin for $SAMPLE"
