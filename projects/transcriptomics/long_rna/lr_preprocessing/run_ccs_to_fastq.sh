#!/usr/bin/env bash
set -euo pipefail

usage(){ echo "Usage: $0 -i <ccs.bam> -o <output.fastq.gz>"; exit 1; }

while getopts "i:o:" opt; do
  case ${opt} in
    i) CCS=$OPTARG ;;
    o) OUTPUT=$OPTARG ;;
    *) usage ;;
  esac
done
[[ -z "$CCS" || -z "$OUTPUT" ]] && usage

mkdir -p "$(dirname "$OUTPUT")"

samtools fastq "$CCS" | gzip > "$OUTPUT"

echo "✅ FASTQ export → $OUTPUT"
