#!/usr/bin/env bash
set -euo pipefail

usage() { echo "Usage: $0 -i <input.fastq.gz> -o <output_dir>"; exit 1; }

while getopts "i:o:" opt; do
  case ${opt} in
    i) INPUT=$OPTARG ;;
    o) OUTDIR=$OPTARG ;;
    *) usage ;;
  esac
done

[[ -z "$INPUT" || -z "$OUTDIR" ]] && usage
mkdir -p "$OUTDIR"

echo "📊 Generating NanoPlot QC..."
NanoPlot --fastq "$INPUT" --loglength --plots hex dot -o "$OUTDIR"

echo "✅ QC report → $OUTDIR"
