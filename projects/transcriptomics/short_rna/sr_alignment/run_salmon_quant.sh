#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 -i <salmon_index> -s <sample_sheet.tsv> -o <output_dir> [-l <libtype>] [-t <threads>]"
  exit 1
}

THREADS=8
LIBTYPE=A
while getopts "i:s:o:l:t:" opt; do
  case $opt in
    i) INDEX=$OPTARG ;;
    s) SHEET=$OPTARG ;;
    o) OUTDIR=$OPTARG ;;
    l) LIBTYPE=$OPTARG ;;
    t) THREADS=$OPTARG ;;
    *) usage ;;
  esac
done

[[ -z "${INDEX:-}" || -z "${SHEET:-}" || -z "${OUTDIR:-}" ]] && usage
[[ ! -d "$INDEX" ]] && { echo "❌ Salmon index not found"; exit 1; }

mkdir -p "$OUTDIR"

while IFS=$'\t' read -r ID R1 R2; do
  echo "🎯 Quantifying sample: $ID"
  salmon quant -i "$INDEX" -l "$LIBTYPE" -1 "$R1" -2 "$R2" -p "$THREADS" -o "$OUTDIR/${ID}_quant"
done < <(tail -n +2 "$SHEET")

echo "✅ Salmon quantification complete. Outputs in $OUTDIR"
