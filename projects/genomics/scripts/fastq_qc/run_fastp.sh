#!/bin/bash
# fastp adapter trimming and QC
# Input: paired-end FASTQ files
# Output: trimmed FASTQ files + HTML report

FQ1=$1
FQ2=$2
OUT1=$3
OUT2=$4
REPORT=$5

fastp -i $FQ1 -I $FQ2 -o $OUT1 -O $OUT2 -h $REPORT -q 20 -l 36
