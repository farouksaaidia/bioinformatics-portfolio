#!/bin/bash
# align_hisat2.sh → Align single bulk RNA-seq sample with HISAT2
# Input: FASTQ R1/R2, reference genome prefix, output folder
# Output: Sorted BAM alignment

FASTQ1=$1
FASTQ2=$2
REF_PREFIX=$3
OUT_DIR=$4

mkdir -p $OUT_DIR

hisat2 -x $REF_PREFIX -1 $FASTQ1 -2 $FASTQ2 -p 8 | samtools sort -o $OUT_DIR/sample_hisat2.bam
samtools index $OUT_DIR/sample_hisat2.bam
