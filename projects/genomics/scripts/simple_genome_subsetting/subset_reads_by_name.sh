#!/bin/bash
# subset_reads_by_name.sh → Extract reads from BAM/FASTQ by read IDs
# Usage: ./subset_reads_by_name.sh input.bam read_ids.txt out.bam
# Input: BAM/FASTQ file + text file with read IDs
# Output: BAM/FASTQ with only those reads

IN=$1
READLIST=$2
OUT=$3

# Detect if input is BAM or FASTQ
if [[ $IN == *.bam ]]; then
    samtools view -N $READLIST -b $IN -o $OUT
else
    seqtk subseq $IN $READLIST > $OUT
fi
