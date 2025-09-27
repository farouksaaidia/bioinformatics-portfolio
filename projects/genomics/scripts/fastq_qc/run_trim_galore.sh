#!/bin/bash
# Trim Galore adapter and quality trimming
# Input: paired-end FASTQ files
# Output: trimmed FASTQ files

FQ1=$1
FQ2=$2
OUT_DIR=$3

trim_galore --paired --quality 20 --fastqc -o $OUT_DIR $FQ1 $FQ2
