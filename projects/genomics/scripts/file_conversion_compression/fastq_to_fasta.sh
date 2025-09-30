#!/usr/bin/env bash
# Convert FASTQ to FASTA
# Usage: ./fastq_to_fasta.sh input.fastq.gz output.fasta

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input.fastq(.gz) output.fasta"
    exit 1
fi

infile=$1
outfile=$2

# zcat if gzipped, else cat
if [[ $infile == *.gz ]]; then
    zcat "$infile" | awk 'NR%4==1{print ">"substr($0,2)} NR%4==2{print}' > "$outfile"
else
    awk 'NR%4==1{print ">"substr($0,2)} NR%4==2{print}' "$infile" > "$outfile"
fi

echo "✅ FASTQ converted to FASTA: $outfile"
