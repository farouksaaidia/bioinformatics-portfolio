#!/bin/bash
# Run RSeQC checks on BAM file(s)

INPUT=$1
REF_GTF=$2
OUTPUT=${3:-rseqc_results}

mkdir -p "$OUTPUT"

if [ -d "$INPUT" ]; then
    for bam in "$INPUT"/*.bam; do
        prefix=$(basename "$bam" .bam)
        junction_saturation.py -i "$bam" -r "$REF_GTF" -o "$OUTPUT/$prefix"
        geneBody_coverage.py -i "$bam" -r "$REF_GTF" -o "$OUTPUT/$prefix"
    done
else
    prefix=$(basename "$INPUT" .bam)
    junction_saturation.py -i "$INPUT" -r "$REF_GTF" -o "$OUTPUT/$prefix"
    geneBody_coverage.py -i "$INPUT" -r "$REF_GTF" -o "$OUTPUT/$prefix"
fi
