#!/usr/bin/env bash
set -euo pipefail

DB=$1
BAM=$2
GTF=$3
OUTDIR=$4

mkdir -p "$OUTDIR"

talon --f "$BAM" --db "$DB" --t 8 --build hg38 --o "$OUTDIR/talon"
talon_abundance --db "$DB" --build hg38 --o "$OUTDIR/talon_abundance"

echo "📚 TALON annotation & quantification complete → $OUTDIR"
