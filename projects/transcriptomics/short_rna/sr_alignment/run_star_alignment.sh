#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 -r <reference_dir> -s <sample_sheet.tsv> -o <output_dir> [-t <threads>]"
  exit 1
}

THREADS=8
while getopts "r:s:o:t:" opt; do
  case $opt in
    r) REF=$OPTARG ;;
    s) SHEET=$OPTARG ;;
    o) OUTDIR=$OPTARG ;;
    t) THREADS=$OPTARG ;;
    *) usage ;;
  esac
done

[[ -z "${REF:-}" || -z "${SHEET:-}" || -z "${OUTDIR:-}" ]] && usage
[[ ! -d "$REF" ]] && { echo "❌ Reference directory not found"; exit 1; }
[[ ! -f "$SHEET" ]] && { echo "❌ Sample sheet not found"; exit 1; }

mkdir -p "$OUTDIR"

while IFS=$'\t' read -r ID R1 R2; do
  echo "🚀 Aligning sample: $ID"
  STAR --runThreadN "$THREADS" \
       --genomeDir "$REF" \
       --readFilesIn "$R1" "$R2" \
       --readFilesCommand zcat \
       --outFileNamePrefix "$OUTDIR/${ID}_" \
       --outSAMtype BAM SortedByCoordinate
done < <(tail -n +2 "$SHEET")

echo "✅ STAR alignment complete. BAMs saved in $OUTDIR"
