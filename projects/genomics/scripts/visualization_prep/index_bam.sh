#!/bin/bash
# index_bam.sh → Create BAM index for visualization
# Input: sorted BAM
# Output: BAM index (.bai)

BAM=$1

samtools index $BAM
