#!/bin/bash
# BLASTN workflow for DNA reads (single-end mock usage)

REF=~/bioinformatics-portfolio/projects/genomics/data/reference/human_g1k_v37.fasta
READS=~/bioinformatics-portfolio/projects/genomics/data/real_trial/SRR1910373_1.fastq.gz
OUTDIR=~/bioinformatics-portfolio/projects/genomics/results/real_results
mkdir -p $OUTDIR

# Prepare BLAST database (only once needed)
makeblastdb -in $REF -dbtype nucl -out $REF.db

# Align small reads (convert FASTQ → FASTA first for BLAST)
zcat $READS | seqtk seq -A - > $OUTDIR/reads.fasta
blastn -query $OUTDIR/reads.fasta -db $REF.db -out $OUTDIR/alignment_blast.txt -outfmt 6

# BLAST doesn't produce BAM, but we store the tabular results
