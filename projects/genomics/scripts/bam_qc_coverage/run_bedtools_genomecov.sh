#!/usr/bin/env bash
# Run bedtools genomecov for genome-wide coverage profile

if [ $# -lt 2 ]; then
    echo "Usage: $0 <input.bam> <reference.fai>"
    exit 1
fi

BAM=$1
FAI=$2
OUT=${BAM%.bam}_genomecov.bedgraph

bedtools genomecov -ibam "$BAM" -g "$FAI" -bg > "$OUT"

echo "Genome coverage written to $OUT"
