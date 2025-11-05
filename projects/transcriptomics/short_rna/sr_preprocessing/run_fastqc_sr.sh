#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 -i <input_fastq_dir> -o <output_dir> [-t <threads>]"
  exit 1
}

THREADS=4
while getopts "i:o:t:" opt; do
  case $opt in
    i) INDIR=$OPTARG ;;
    o) OUTDIR=$OPTARG ;;
    t) THREADS=$OPTARG ;;
    *) usage ;;
  esac
done

[[ -z "${INDIR:-}" || -z "${OUTDIR:-}" ]] && usage
mkdir -p "$OUTDIR"

echo "🔍 Running FastQC on all FASTQ files in $INDIR"
fastqc -t "$THREADS" -o "$OUTDIR" "$INDIR"/*.fastq* || { echo "❌ FastQC failed"; exit 1; }
echo "✅ FastQC complete. Results in $OUTDIR"
