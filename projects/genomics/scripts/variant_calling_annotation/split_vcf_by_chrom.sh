#!/bin/bash
# Split large VCF into per-chromosome files
VCF=$1
OUTDIR=$2

if [[ ! -f $VCF ]]; then
    echo "Error: $VCF not found."
    exit 1
fi

mkdir -p $OUTDIR

bcftools index $VCF
for CHR in $(bcftools query -l $VCF | cut -f1 | sort | uniq); do
    bcftools view -r $CHR $VCF -O v -o $OUTDIR/${CHR}.vcf
    echo "Saved $OUTDIR/${CHR}.vcf"
done
