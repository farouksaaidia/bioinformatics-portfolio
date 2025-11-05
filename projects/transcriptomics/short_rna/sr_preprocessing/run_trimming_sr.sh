#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 -i <input_dir> -o <output_dir> -a <adapters.fa> -m <mode:SE|PE> [-t <threads>]"
  exit 1
}

THREADS=4
while getopts "i:o:a:m:t:" opt; do
  case $opt in
    i) INDIR=$OPTARG ;;
    o) OUTDIR=$OPTARG ;;
    a) ADAPTERS=$OPTARG ;;
    m) MODE=$OPTARG ;;
    t) THREADS=$OPTARG ;;
    *) usage ;;
  esac
done

[[ -z "${INDIR:-}" || -z "${OUTDIR:-}" || -z "${ADAPTERS:-}" || -z "${MODE:-}" ]] && usage
mkdir -p "$OUTDIR"

if [[ "$MODE" == "SE" ]]; then
  echo "✂️  Running Cutadapt (single-end)"
  for f in "$INDIR"/*fastq*; do
    base=$(basename "$f")
    cutadapt -a file:"$ADAPTERS" -q 20 -j "$THREADS" -o "$OUTDIR/$base" "$f"
  done
elif [[ "$MODE" == "PE" ]]; then
  echo "✂️  Running Cutadapt (paired-end)"
  for f in "$INDIR"/*_R1.fastq*; do
    base=$(basename "$f" _R1.fastq.gz)
    cutadapt -a file:"$ADAPTERS" -A file:"$ADAPTERS" -q 20 -j "$THREADS" \
      -o "$OUTDIR/${base}_R1.trimmed.fastq.gz" -p "$OUTDIR/${base}_R2.trimmed.fastq.gz" \
      "$INDIR/${base}_R1.fastq.gz" "$INDIR/${base}_R2.fastq.gz"
  done
else
  echo "❌ Mode must be SE or PE"; exit 1
fi

echo "✅ Trimming complete. Results in $OUTDIR"
