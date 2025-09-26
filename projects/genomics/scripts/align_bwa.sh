#!/bin/bash
# BWA alignment workflow: FASTQ → SAM → BAM → Sorted BAM → Index

# Input files
REF=~/bioinformatics-portfolio/projects/genomics/data/reference/human_g1k_v37.fasta
READS1=~/bioinformatics-portfolio/projects/genomics/data/real_trial/SRR1910373_1.fastq.gz
READS2=~/bioinformatics-portfolio/projects/genomics/data/real_trial/SRR1910373_2.fastq.gz

# Output directory
OUTDIR=~/bioinformatics-portfolio/projects/genomics/results/real_results
mkdir -p $OUTDIR

# Run BWA MEM
bwa mem $REF $READS1 $READS2 > $OUTDIR/alignment_bwa.sam

# SAM → BAM
samtools view -bS $OUTDIR/alignment_bwa.sam > $OUTDIR/alignment_bwa.bam

# Sort BAM
samtools sort $OUTDIR/alignment_bwa.bam -o $OUTDIR/alignment_bwa.sorted.bam

# Index BAM
samtools index $OUTDIR/alignment_bwa.sorted.bam

# Cleanup big intermediate
rm $OUTDIR/alignment_bwa.sam $OUTDIR/alignment_bwa.bam
