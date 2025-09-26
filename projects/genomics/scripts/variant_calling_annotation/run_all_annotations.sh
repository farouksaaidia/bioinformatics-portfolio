#!/bin/bash
# run_all_annotations.sh → Annotate a VCF with all annotation scripts
# Input: VCF file, output directory
# Output: fully annotated VCF

VCF=$1
OUTDIR=$2

mkdir -p $OUTDIR

# Sequentially annotate
./annotate_vep.sh $VCF $OUTDIR/annotated_vep.vcf
./annotate_annovar.sh $OUTDIR/annotated_vep.vcf $OUTDIR/annotated_annovar.vcf
./annotate_variants_snpeff.sh $OUTDIR/annotated_annovar.vcf $OUTDIR/annotated_snpeff.vcf
./annotate_bcftools_csq.sh $OUTDIR/annotated_snpeff.vcf $OUTDIR/annotated_bcftools.vcf
./annotate_population_frequency.sh $OUTDIR/annotated_bcftools.vcf $OUTDIR/annotated_population.vcf
./extract_clinical_variants.sh $OUTDIR/annotated_population.vcf $OUTDIR/annotated_clinical.vcf
./annotate_conservation_scores.sh $OUTDIR/annotated_clinical.vcf $OUTDIR/annotated_final.vcf

echo "All annotations complete. Final annotated VCF: $OUTDIR/annotated_final.vcf"
