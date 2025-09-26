#!/bin/bash
# Annotate VCF using Ensembl VEP
# Input: VCF, reference genome, cache
# Output: annotated VCF

VCF=$1
OUT=$2
CACHE_DIR=${3:-$HOME/.vep/cache}

vep -i $VCF \
    --cache \
    --dir_cache $CACHE_DIR \
    --vcf \
    --force_overwrite \
    -o $OUT \
    --everything
