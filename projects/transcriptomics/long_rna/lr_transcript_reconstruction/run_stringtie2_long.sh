#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 -b input.bam -a annotation.gtf -o out_dir -t threads"
  exit 1
}

THREADS=8
while getopts "b:a:o:t:" opt; do
  case $opt in
    b) BAM=$OPTARG ;;
    a) GTF=$OPTARG ;;
    o) OUTDIR=$OPTARG ;;
    t) THREADS=$OPTARG ;;
    *) usage ;;
  esac
done

mkdir -p "$OUTDIR"

stringtie "$BAM" -L -G "$GTF" -o "$OUTDIR/stringtie.gtf" -A "$OUTDIR/gene_abundance.tsv" -p "$THREADS"

echo "📊 StringTie2 quantification complete → $OUTDIR"
