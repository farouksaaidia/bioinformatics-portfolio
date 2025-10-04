#!/bin/bash
# =========================================
# run_starsolo.sh
# Preprocessing pipeline for scRNA-seq using STARsolo
# =========================================

set -euo pipefail

usage() {
    echo "Usage: $0 -i <sample_id> -1 <R1.fastq.gz> -2 <R2.fastq.gz> -r <reference_dir> -o <output_dir>"
    exit 1
}

while getopts ":i:1:2:r:o:h" opt; do
    case $opt in
        i) SAMPLE=$OPTARG ;;
        1) R1=$OPTARG ;;
        2) R2=$OPTARG ;;
        r) REF=$OPTARG ;;
        o) OUTDIR=$OPTARG ;;
        h) usage ;;
    esac
done

if [[ -z "${SAMPLE:-}" || -z "${R1:-}" || -z "${R2:-}" || -z "${REF:-}" || -z "${OUTDIR:-}" ]]; then
    echo "[ERROR] Missing arguments"; usage
fi

echo "[$(date)] Running STARsolo for $SAMPLE"
STAR --runThreadN 8 --genomeDir "$REF" --readFilesIn "$R1" "$R2" --readFilesCommand zcat --soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 --outFileNamePrefix "$OUTDIR/$SAMPLE."
echo "[$(date)] Finished STARsolo for $SAMPLE"
