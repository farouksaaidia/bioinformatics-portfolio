#!/usr/bin/env bash
set -euo pipefail

usage() { echo "Usage: $0 -i <subreads.bam> -o <ccs.bam>"; exit 1; }

while getopts "i:o:" opt; do
  case ${opt} in
    i) SUBREADS=$OPTARG ;;
    o) CCS=$OPTARG ;;
    *) usage ;;
  esac
done
[[ -z "$SUBREADS" || -z "$CCS" ]] && usage
mkdir -p "$(dirname "$CCS")"

echo "🔄 Creating CCS HiFi reads..."
ccs "$SUBREADS" "$CCS" --numThreads=8

echo "✅ CCS reads → $CCS"
