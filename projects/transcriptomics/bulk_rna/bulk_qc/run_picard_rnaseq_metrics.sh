#!/bin/bash
# Run Picard CollectRnaSeqMetrics on BAM file(s)

INPUT=$1
REF_FLAT=$2   # RefFlat file
RIBO_INTERVALS=$3 # Ribosomal intervals list
OUTPUT=${4:-picard_rnaseq_metrics}

mkdir -p "$OUTPUT"

if [ -d "$INPUT" ]; then
    for bam in "$INPUT"/*.bam; do
        prefix=$(basename "$bam" .bam)
        picard CollectRnaSeqMetrics \
            I="$bam" \
            O="$OUTPUT/${prefix}_rnaseq_metrics.txt" \
            REF_FLAT="$REF_FLAT" \
            RIBOSOMAL_INTERVALS="$RIBO_INTERVALS"
    done
else
    prefix=$(basename "$INPUT" .bam)
    picard CollectRnaSeqMetrics \
        I="$INPUT" \
        O="$OUTPUT/${prefix}_rnaseq_metrics.txt" \
        REF_FLAT="$REF_FLAT" \
        RIBOSOMAL_INTERVALS="$RIBO_INTERVALS"
fi
