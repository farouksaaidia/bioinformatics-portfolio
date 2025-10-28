#!/usr/bin/env bash
set -euo pipefail

INPUT_LOG=$1

grep -E "mapped|primary" "$INPUT_LOG" > alignment_summary.txt

echo "📊 Extracted key mapping statistics → alignment_summary.txt"
