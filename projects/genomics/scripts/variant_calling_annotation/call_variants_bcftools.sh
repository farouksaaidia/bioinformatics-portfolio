#!/bin/bash
# BCFtools variant calling
# Input: sorted BAM, reference FASTA
# Output: raw VCF

BAM=$1
REF=$2
OUT=$3
THREADS=${4:-4}

# Index the reference if not already indexed
samtools faidx $REF

# Generate mpileup and call variants
bcftools mpileup -f $REF -Ou -b <(echo $BAM) -a FORMAT/DP | \
    bcftools call -mv -Ov -o $OUT
