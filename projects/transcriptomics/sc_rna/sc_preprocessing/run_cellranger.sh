#!/bin/bash
# =========================================
# run_cellranger.sh
# Preprocessing pipeline for scRNA-seq using CellRanger
# Supports single or multiple samples
# =========================================

set -euo pipefail

usage() {
    echo "Usage: $0 -i <sample_id> -f <fastq_path> -r <reference_path> -o <output_dir>"
    echo "Example (single sample):"
    echo "  $0 -i Sample1 -f ./fastq -r ./refdata-gex-GRCh38-2020-A -o ./cellranger_out"
    echo "Example (multiple samples via manifest.csv):"
    echo "  while IFS=, read SAMPLE FASTQ REF OUTDIR; do"
    echo "      $0 -i \$SAMPLE -f \$FASTQ -r \$REF -o \$OUTDIR"
    echo "  done < manifest.csv"
    exit 1
}

while getopts ":i:f:r:o:h" opt; do
    case $opt in
        i) SAMPLE=$OPTARG ;;
        f) FASTQ=$OPTARG ;;
        r) REF=$OPTARG ;;
        o) OUTDIR=$OPTARG ;;
        h) usage ;;
        \?) echo "Invalid option -$OPTARG"; usage ;;
    esac
done

if [[ -z "${SAMPLE:-}" || -z "${FASTQ:-}" || -z "${REF:-}" || -z "${OUTDIR:-}" ]]; then
    echo "[ERROR] Missing arguments"; usage
fi

echo "[$(date)] Starting CellRanger for sample: $SAMPLE"
cellranger count --id="$SAMPLE" --fastqs="$FASTQ" --transcriptome="$REF" --sample="$SAMPLE" --localcores=8 --localmem=64 --output-dir="$OUTDIR"
echo "[$(date)] Finished CellRanger for sample: $SAMPLE"
