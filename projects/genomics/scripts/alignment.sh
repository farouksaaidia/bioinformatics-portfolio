#!/bin/bash
# alignment.sh - Align paired-end FASTQ reads to a reference genome using BWA + SAMtools
# Usage: ./alignment.sh reads_1.fastq reads_2.fastq reference.fasta output_prefix

# Input arguments
READ1=$1
READ2=$2
REF=$3
OUT_PREFIX=$4

# Step 1: Index the reference genome
bwa index $REF

# Step 2: Align reads to reference
bwa mem $REF $READ1 $READ2 > ${OUT_PREFIX}.sam

# Step 3: Convert SAM -> BAM
samtools view -S -b ${OUT_PREFIX}.sam > ${OUT_PREFIX}.bam

# Step 4: Sort BAM
samtools sort ${OUT_PREFIX}.bam -o ${OUT_PREFIX}_sorted.bam

# Step 5: Index BAM
samtools index ${OUT_PREFIX}_sorted.bam

echo "Alignment complete! Results:"
echo "- SAM file: ${OUT_PREFIX}.sam"
echo "- BAM file: ${OUT_PREFIX}.bam"
echo "- Sorted BAM: ${OUT_PREFIX}_sorted.bam"
echo "- BAM index: ${OUT_PREFIX}_sorted.bam.bai"

