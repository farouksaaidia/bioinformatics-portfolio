#!/bin/bash
# =========================================
# run_zumis.sh
# Preprocessing pipeline for scRNA-seq using zUMIs
# =========================================

set -euo pipefail

usage() {
    echo "Usage: $0 -1 <R1.fastq.gz> -2 <R2.fastq.gz> -g <genome.fa> -a <annotation.gtf> -o <output_dir>"
    exit 1
}

while getopts ":1:2:g:a:o:h" opt; do
    case $opt in
        1) R1=$OPTARG ;;
        2) R2=$OPTARG ;;
        g) GENOME=$OPTARG ;;
        a) ANNO=$OPTARG ;;
        o) OUTDIR=$OPTARG ;;
        h) usage ;;
    esac
done

if [[ -z "${R1:-}" || -z "${R2:-}" || -z "${GENOME:-}" || -z "${ANNO:-}" || -z "${OUTDIR:-}" ]]; then
    echo "[ERROR] Missing arguments"; usage
fi

echo "[$(date)] Running zUMIs"
zUMIs -c config.yaml -1 "$R1" -2 "$R2" -g "$GENOME" -a "$ANNO" -o "$OUTDIR"
echo "[$(date)] Finished zUMIs"
