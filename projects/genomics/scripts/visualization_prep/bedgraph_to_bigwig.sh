#!/bin/bash
# bedgraph_to_bigwig.sh → Convert bedGraph to BigWig
# Input: bedGraph, genome.chrom.sizes
# Output: BigWig

BEDGRAPH=$1
CHROMSIZES=$2
OUT=$3

bedGraphToBigWig $BEDGRAPH $CHROMSIZES $OUT
