#!/bin/bash
# subset_bam_region.sh → Extract reads from BAM in a region
# Input: BAM, region (chr:start-end)
# Output: subset BAM

BAM=$1
REGION=$2
OUT=$3

samtools view -b $BAM $REGION -o $OUT
