#!/bin/bash
# run_all_callers.sh → Call variants with GATK, bcftools, FreeBayes and merge them
# Input: sorted BAM, reference FASTA, output directory
# Output: separate VCFs + merged VCF

BAM=$1
REF=$2
OUTDIR=$3

mkdir -p $OUTDIR

# Call variants
./call_variants_gatk.sh $BAM $REF $OUTDIR/gatk.vcf
./call_variants_bcftools.sh $BAM $REF $OUTDIR/bcftools.vcf
./call_freebayes.sh $BAM $REF $OUTDIR/freebayes.vcf

# Merge VCFs into a single non-redundant VCF
bcftools merge -m none $OUTDIR/gatk.vcf $OUTDIR/bcftools.vcf $OUTDIR/freebayes.vcf -Oz -o $OUTDIR/merged.vcf.gz
bcftools index $OUTDIR/merged.vcf.gz

echo "Variant calling complete. Separate VCFs and merged VCF generated in $OUTDIR"
