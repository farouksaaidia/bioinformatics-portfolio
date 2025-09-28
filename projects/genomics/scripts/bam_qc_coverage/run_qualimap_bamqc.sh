#!/usr/bin/env bash
# Run Qualimap BAM QC

if [ $# -lt 2 ]; then
    echo "Usage: $0 <input.bam> <reference.fa>"
    exit 1
fi

BAM=$1
REF=$2
OUTDIR=${BAM%.bam}_qualimap

qualimap bamqc -bam "$BAM" -gff "$REF" -outdir "$OUTDIR" -outformat HTML,PDF

echo "Qualimap BAM QC written to $OUTDIR"
