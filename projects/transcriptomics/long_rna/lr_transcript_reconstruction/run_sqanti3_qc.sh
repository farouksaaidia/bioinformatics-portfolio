#!/usr/bin/env bash
set -euo pipefail

GTF=$1
GENOME=$2
ANN=$3
OUTDIR=$4

mkdir -p "$OUTDIR"

sqanti3_qc.py "$GTF" "$ANN" "$GENOME" --output_dir "$OUTDIR" --fl_count collapsed_abundance.tsv

echo "✅ SQANTI3 QC and isoform classification done → $OUTDIR"
