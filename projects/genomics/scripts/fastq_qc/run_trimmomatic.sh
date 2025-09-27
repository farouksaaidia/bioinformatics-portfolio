#!/bin/bash
# Trimmomatic adapter and quality trimming
# Input: paired-end FASTQ files
# Output: trimmed FASTQ files

FQ1=$1
FQ2=$2
OUT1P=$3
OUT1U=$4
OUT2P=$5
OUT2U=$6

trimmomatic PE -threads 4 $FQ1 $FQ2 $OUT1P $OUT1U $OUT2P $OUT2U ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
