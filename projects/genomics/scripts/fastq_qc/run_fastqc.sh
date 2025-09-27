#!/bin/bash
# run_fastqc.sh → Run FastQC on raw FASTQ files
# Usage: ./run_fastqc.sh file1.fastq.gz [file2.fastq.gz ...]

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <fastq_file1> [fastq_file2 ...]"
    exit 1
fi

mkdir -p fastqc_results

for fq in "$@"; do
    echo "Running FastQC on $fq..."
    fastqc -o fastqc_results -t 4 "$fq"
done

echo "FastQC analysis complete. Results in fastqc_results/"
