#!/bin/bash
# Annotate VCF using ANNOVAR
# Input: VCF, ANNOVAR database
# Output: annotated table

VCF=$1
DB=$2
OUT=$3

# Convert VCF to ANNOVAR input
convert2annovar.pl -format vcf4 $VCF -outfile ${OUT%.txt}.avinput

# Annotate variants
table_annovar.pl ${OUT%.txt}.avinput $DB \
    -buildver hg38 \
    -out $OUT \
    -remove \
    -protocol refGene,cytoBand,dbnsfp41a,clinvar_20220320 \
    -operation g,r,f,f \
    -nastring . \
    -vcfinput
