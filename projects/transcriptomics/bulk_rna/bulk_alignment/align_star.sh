#!/bin/bash
# align_star.sh → Align single bulk RNA-seq sample with STAR
# Input: FASTQ R1/R2, reference genome dir, output folder
# Output: Sorted BAM alignment

FASTQ1=$1
FASTQ2=$2
REF_DIR=$3
OUT_DIR=$4

mkdir -p $OUT_DIR

STAR --genomeDir $REF_DIR \
     --readFilesIn $FASTQ1 $FASTQ2 \
     --readFilesCommand zcat \
     --runThreadN 8 \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix $OUT_DIR/sample_
