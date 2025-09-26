#!/bin/bash
# Bowtie2 alignment workflow

REF=~/bioinformatics-portfolio/projects/genomics/data/reference/human_g1k_v37
READS1=~/bioinformatics-portfolio/projects/genomics/data/real_trial/SRR1910373_1.fastq.gz
READS2=~/bioinformatics-portfolio/projects/genomics/data/real_trial/SRR1910373_2.fastq.gz
OUTDIR=~/bioinformatics-portfolio/projects/genomics/results/real_results
mkdir -p $OUTDIR

bowtie2 -x $REF -1 $READS1 -2 $READS2 -S $OUTDIR/alignment_bt2.sam

samtools view -bS $OUTDIR/alignment_bt2.sam > $OUTDIR/alignment_bt2.bam
samtools sort $OUTDIR/alignment_bt2.bam -o $OUTDIR/alignment_bt2.sorted.bam
samtools index $OUTDIR/alignment_bt2.sorted.bam

rm $OUTDIR/alignment_bt2.sam $OUTDIR/alignment_bt2.bam
