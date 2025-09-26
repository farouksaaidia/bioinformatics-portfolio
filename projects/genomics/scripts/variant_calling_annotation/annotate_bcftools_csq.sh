#!/bin/bash
# Annotate VCF using bcftools csq
# Input: VCF, reference FASTA, gene annotation (GFF/GTF)
# Output: annotated VCF

VCF=$1
REF=$2
GTF=$3
OUT=$4

# Index reference if needed
samtools faidx $REF

# Annotate variants
bcftools csq -f $REF -g $GTF $VCF -o $OUT -Ob
bcftools view $OUT -Ov -o ${OUT%.bcf}.vcf
