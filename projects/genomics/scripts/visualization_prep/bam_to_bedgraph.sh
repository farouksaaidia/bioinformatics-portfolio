#!/bin/bash
# bam_to_bedgraph.sh → Convert BAM to bedGraph for coverage visualization
# Input: sorted BAM
# Output: bedGraph

BAM=$1
OUT=$2

genomeCoverageBed -ibam $BAM -bg > $OUT
