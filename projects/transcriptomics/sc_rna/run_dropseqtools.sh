#!/bin/bash
# =========================================
# run_dropseqtools.sh
# Preprocessing pipeline for scRNA-seq using Drop-seq Tools
# =========================================

set -euo pipefail

usage() {
    echo "Usage: $0 -r <R1.fastq.gz> -s <R2.fastq.gz> -g <genome.fa> -a <annotation.gtf> -o <output_dir>"
    exit 1
}

while getopts ":r:s:g:a:o:h" opt; do
    case $opt in
        r) R1=$OPTARG ;;
        s) R2=$OPTARG ;;
        g) GENOME=$OPTARG ;;
        a) ANNO=$OPTARG ;;
        o) OUTDIR=$OPTARG ;;
        h) usage ;;
    esac
done

if [[ -z "${R1:-}" || -z "${R2:-}" || -z "${GENOME:-}" || -z "${ANNO:-}" || -z "${OUTDIR:-}" ]]; then
    echo "[ERROR] Missing arguments"; usage
fi

echo "[$(date)] Running Drop-seq Tools pipeline"
mkdir -p "$OUTDIR"
dropseq_tools TagBamWithReadSequenceExtended I=unaligned.bam O=cell_tagged.bam SUMMARY=summary.txt BASE_RANGE=1-12 BASE_QUALITY=10 TAG_NAME=XC
echo "[$(date)] Finished Drop-seq Tools"
