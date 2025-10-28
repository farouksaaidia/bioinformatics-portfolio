#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 -r <reference.fasta> -i <input.fastq.gz or list.txt> -o <output_dir> [-t threads]"
  exit 1
}

THREADS=8
while getopts "r:i:o:t:" opt; do
  case $opt in
    r) REF=$OPTARG ;;
    i) INPUT=$OPTARG ;;
    o) OUTDIR=$OPTARG ;;
    t) THREADS=$OPTARG ;;
    *) usage ;;
  esac
done

[[ -z "${REF:-}" || -z "${INPUT:-}" || -z "${OUTDIR:-}" ]] && usage

mkdir -p "$OUTDIR"
[ ! -f "$REF" ] && { echo "Reference not found"; exit 1; }

map_sample() {
  fq=$1
  sample=$(basename "$fq" .fastq.gz)
  bam="$OUTDIR/${sample}.bam"

  minimap2 -ax splice -uf -k14 -t "$THREADS" "$REF" "$fq" \
    | samtools sort -@ "$THREADS" -o "$bam"
  samtools index "$bam"

  echo "✅ Aligned sample: $sample"
}

if [[ -f "$INPUT" && "$INPUT" == *.txt ]]; then
  while read fq; do map_sample "$fq"; done < "$INPUT"
else
  map_sample "$INPUT"
fi

echo "🎯 Alignment complete → $OUTDIR"
