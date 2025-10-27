#!/usr/bin/env bash
set -euo pipefail

usage() { echo "Usage: $0 -i <input.fastq.gz> -o <output.fastq.gz>"; exit 1; }

while getopts "i:o:" opt; do
  case ${opt} in
    i) INPUT=$OPTARG ;;
    o) OUTPUT=$OPTARG ;;
    *) usage ;;
  esac
done

[[ -z "$INPUT" || -z "$OUTPUT" ]] && usage
mkdir -p "$(dirname $OUTPUT)"

echo "✂️  Trimming adapters with Porechop..."
porechop -i "$INPUT" -o "$OUTPUT" --threads 8

echo "✅ Trimmed reads → $OUTPUT"
