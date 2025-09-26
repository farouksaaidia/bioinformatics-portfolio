#!/bin/bash
# GATK HaplotypeCaller variant calling (realistic)
# Input: sorted BAM, reference FASTA, output directory
# Output: raw VCF

BAM=$1
REF=$2
OUT=$3
THREADS=${4:-4}

# Index reference if needed
if [ ! -f ${REF}.fai ]; then
    samtools faidx $REF
fi

# Create sequence dictionary if not exists
DICT=${REF%.fa}.dict
if [ ! -f $DICT ]; then
    gatk CreateSequenceDictionary -R $REF -O $DICT
fi

# Call variants
gatk HaplotypeCaller \
    -R $REF \
    -I $BAM \
    -O $OUT \
    --native-pair-hmm-threads $THREADS
