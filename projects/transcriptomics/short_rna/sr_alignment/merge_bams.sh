#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 -i <bam_dir> -o <merged_bam>"
  exit 1
}

while getopts "i:o:" opt; do
  case $opt in
    i) BAMDIR=$OPTARG ;;
    o) OUTPUT=$OPTARG ;;
    *) usage ;;
  esac
done

[[ -z "${BAMDIR:-}" || -z "${OUTPUT:-}" ]] && usage
[[ ! -d "$BAMDIR" ]] && { echo "❌ BAM directory not found"; exit 1; }

samtools merge -f "$OUTPUT" "$BAMDIR"/*.bam
samtools index "$OUTPUT"

echo "✅ Merged BAM created at $OUTPUT"
