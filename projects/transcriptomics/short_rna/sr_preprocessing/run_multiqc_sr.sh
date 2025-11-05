#!/usr/bin/env bash
set -euo pipefail

usage() { echo "Usage: $0 -i <fastqc_results_dir> -o <output_dir>"; exit 1; }

while getopts "i:o:" opt; do
  case $opt in
    i) INDIR=$OPTARG ;;
    o) OUTDIR=$OPTARG ;;
    *) usage ;;
  esac
done

[[ -z "${INDIR:-}" || -z "${OUTDIR:-}" ]] && usage
mkdir -p "$OUTDIR"

echo "📊 Generating MultiQC summary..."
multiqc "$INDIR" -o "$OUTDIR" || { echo "❌ MultiQC failed"; exit 1; }
echo "✅ MultiQC report saved to $OUTDIR"
