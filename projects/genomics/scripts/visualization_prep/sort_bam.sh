#!/bin/bash
# sort_bam.sh → Sort BAM file by coordinate
# Input: unsorted BAM
# Output: sorted BAM

BAM=$1
OUT=$2

samtools sort -o $OUT $BAM
