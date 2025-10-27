#!/usr/bin/env bash
set -euo pipefail

usage() { echo "Usage: $0 -i <fast5_dir> -o <output_dir> -c <config>"; exit 1; }

CONFIG="rna_r9.4.1.cfg"

while getopts "i:o:c:" opt; do
  case ${opt} in
    i) INPUT=$OPTARG ;;
    o) OUTPUT=$OPTARG ;;
    c) CONFIG=$OPTARG ;;
    *) usage ;;
  esac
done

[[ -z "$INPUT" || -z "$OUTPUT" ]] && usage
mkdir -p "$OUTPUT"

echo "🔄 Running Guppy basecalling..."
guppy_basecaller \
  -i "$INPUT" \
  -s "$OUTPUT" \
  -c "$CONFIG" \
  --recursive \
  --compress_fastq \
  --num_callers 4 \
  --gpu_runners_per_device 4 \
  --device cuda:all || { echo "❌ Guppy failed"; exit 1; }

echo "✅ Basecalling done → $OUTPUT"
