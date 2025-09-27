#!/bin/bash
# Cutadapt adapter trimming
# Input: paired-end FASTQ files
# Output: trimmed FASTQ files

FQ1=$1
FQ2=$2
OUT1=$3
OUT2=$4

cutadapt -q 20 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o $OUT1 -p $OUT2 $FQ1 $FQ2
