#!/bin/bash
# Run full QC pipeline (FastQC, MultiQC, RSeQC, Picard, Preseq)

FASTQ_DIR=$1
BAM_DIR=$2
REF_GTF=$3
REF_FLAT=$4
RIBO_INTERVALS=$5
OUTDIR=${6:-all_qc_results}

mkdir -p "$OUTDIR"

echo "Running FastQC..."
./run_fastqc.sh "$FASTQ_DIR" "$OUTDIR/fastqc"

echo "Running MultiQC..."
./run_multiqc.sh "$OUTDIR/fastqc" "$OUTDIR/multiqc"

echo "Running RSeQC..."
./run_rseqc.sh "$BAM_DIR" "$REF_GTF" "$OUTDIR/rseqc"

echo "Running Picard CollectRnaSeqMetrics..."
./run_picard_rnaseq_metrics.sh "$BAM_DIR" "$REF_FLAT" "$RIBO_INTERVALS" "$OUTDIR/picard"

echo "Running Preseq..."
./run_preseq.sh "$BAM_DIR" "$OUTDIR/preseq"

echo "All QC completed. Results in $OUTDIR"
