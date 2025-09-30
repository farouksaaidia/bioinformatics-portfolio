#!/bin/bash
# align_bowtie2.sh → Align single bulk RNA-seq sample with Bowtie2
# Input: FASTQ R1/R2, reference genome prefix, output folder
# Output: Sorted BAM alignment

FASTQ1=$1
FASTQ2=$2
REF_PREFIX=$3
OUT_DIR=$4

mkdir -p $OUT_DIR

bowtie2 -x $REF_PREFIX -1 $FASTQ1 -2 $FASTQ2 -p 8 | samtools sort -o $OUT_DIR/sample_bowtie2.bam
samtools index $OUT_DIR/sample_bowtie2.bam
