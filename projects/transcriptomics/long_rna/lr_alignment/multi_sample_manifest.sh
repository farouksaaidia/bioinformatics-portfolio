#!/usr/bin/env bash
set -euo pipefail

MANIFEST=$1
OUTDIR=$2

while IFS=$'\t' read -r fq ref; do
  echo "▶ Running: $fq using $ref"
  ./run_minimap2_align.sh -i "$fq" -r "$ref" -o "$OUTDIR"
done < "$MANIFEST"
