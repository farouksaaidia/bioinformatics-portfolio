#!/bin/bash
# Minimap2 alignment workflow

REF=~/bioinformatics-portfolio/projects/genomics/data/reference/human_g1k_v37.fasta
READS1=~/bioinformatics-portfolio/projects/genomics/data/real_trial/SRR1910373_1.fastq.gz
READS2=~/bioinformatics-portfolio/projects/genomics/data/real_trial/SRR1910373_2.fastq.gz
OUTDIR=~/bioinformatics-portfolio/projects/genomics/results/real_results
mkdir -p $OUTDIR

minimap2 -ax sr $REF $READS1 $READS2 > $OUTDIR/alignment_minimap2.sam

samtools view -bS $OUTDIR/alignment_minimap2.sam > $OUTDIR/alignment_minimap2.bam
samtools sort $OUTDIR/alignment_minimap2.bam -o $OUTDIR/alignment_minimap2.sorted.bam
samtools index $OUTDIR/alignment_minimap2.sorted.bam

rm $OUTDIR/alignment_minimap2.sam $OUTDIR/alignment_minimap2.bam
