#!/bin/bash
# Run RSEM quantification for one or multiple FASTQ files
# Usage:
#   ./run_rsem.sh sample_R1.fastq.gz sample_R2.fastq.gz ref_prefix output_dir
#   ./run_rsem.sh fastq_folder ref_prefix output_dir

set -euo pipefail

INPUT=$1
REF=$2
OUTDIR=$3

mkdir -p "$OUTDIR"

if [ -d "$INPUT" ]; then
    echo "[INFO] Directory detected. Running RSEM on all FASTQ pairs..."
    for R1 in "$INPUT"/*_R1.fastq.gz; do
        R2=${R1/_R1.fastq.gz/_R2.fastq.gz}
        SAMPLE=$(basename "$R1" _R1.fastq.gz)
        rsem-calculate-expression --paired-end -p 8 "$R1" "$R2" "$REF" "$OUTDIR/$SAMPLE"
    done
else
    R1=$INPUT
    R2=${R1/_R1.fastq.gz/_R2.fastq.gz}
    SAMPLE=$(basename "$R1" _R1.fastq.gz)
    rsem-calculate-expression --paired-end -p 8 "$R1" "$R2" "$REF" "$OUTDIR/$SAMPLE"
fi
