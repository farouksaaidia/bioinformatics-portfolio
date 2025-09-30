# Simple Genome Subsetting Scripts

This folder contains helper scripts for extracting subsets of genome or sequencing data.  
These are useful for testing pipelines on smaller datasets, troubleshooting alignments, or performing focused analyses.

## Available Scripts

| Script Name              | Description                                              | Input                         | Output                           | Notes / Tool Options |
|---------------------------|----------------------------------------------------------|-------------------------------|----------------------------------|-----------------------|
| extract_chromosome.sh     | Extract a single chromosome from BAM/FASTA               | BAM/FASTA + chromosome name   | Subset BAM/FASTA                  | Requires samtools     |
| extract_region.sh         | Extract a specific genomic region from BAM/FASTA         | BAM/FASTA + region (chr:start-end) | Subset BAM/FASTA           | Useful for zooming in |
| extract_random_reads.sh   | Extract random reads from FASTQ                          | FASTQ file + number of reads  | Subset FASTQ                      | Useful for testing    |
