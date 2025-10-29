#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 -g genome.fa -a annotation.gtf -b input.bam -o out_dir -t threads"
  exit 1
}

THREADS=8
while getopts "g:a:b:o:t:" opt; do
  case $opt in
    g) GENOME=$OPTARG ;;
    a) GTF=$OPTARG ;;
    b) BAM=$OPTARG ;;
    o) OUTDIR=$OPTARG ;;
    t) THREADS=$OPTARG ;;
    *) usage ;;
  esac
done

[[ -z "${GENOME:-}" || -z "${GTF:-}" || -z "${BAM:-}" || -z "${OUTDIR:-}" ]] && usage

mkdir -p "$OUTDIR"

echo "🧬 Running FLAIR correction..."
flair correct -q "$BAM" -g "$GENOME" -f "$GTF" -o "$OUTDIR/corrected" -t "$THREADS"

echo "🧩 Assembling transcripts..."
flair collapse -g "$GENOME" -r "$OUTDIR/corrected.bam" -f "$GTF" -t "$THREADS" -o "$OUTDIR/isoforms"

echo "🎯 FLAIR pipeline completed → $OUTDIR"
