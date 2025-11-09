#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 -x <hisat2_index_prefix> -s <sample_sheet.tsv> -o <output_dir> [-t <threads>]"
  exit 1
}

THREADS=8
while getopts "x:s:o:t:" opt; do
  case $opt in
    x) INDEX=$OPTARG ;;
    s) SHEET=$OPTARG ;;
    o) OUTDIR=$OPTARG ;;
    t) THREADS=$OPTARG ;;
    *) usage ;;
  esac
done

[[ -z "${INDEX:-}" || -z "${SHEET:-}" || -z "${OUTDIR:-}" ]] && usage
[[ ! -f "${INDEX}.1.ht2" ]] && { echo "❌ HISAT2 index not found"; exit 1; }
[[ ! -f "$SHEET" ]] && { echo "❌ Sample sheet not found"; exit 1; }

mkdir -p "$OUTDIR"

while IFS=$'\t' read -r ID R1 R2; do
  echo "🧬 Running HISAT2 for $ID"
  hisat2 -x "$INDEX" -1 "$R1" -2 "$R2" -p "$THREADS" | samtools sort -o "$OUTDIR/${ID}.bam"
  samtools index "$OUTDIR/${ID}.bam"
done < <(tail -n +2 "$SHEET")

echo "✅ HISAT2 alignment done. BAMs saved in $OUTDIR"
