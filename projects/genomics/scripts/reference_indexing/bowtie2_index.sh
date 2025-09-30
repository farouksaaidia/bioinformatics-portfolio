#!/bin/bash
# Script: bowtie2_index.sh
# Description: Index a reference FASTA for Bowtie2 alignment
# Usage: ./bowtie2_index.sh reference.fasta index_prefix

if [ $# -ne 2 ]; then
    echo "Usage: $0 reference.fasta index_prefix"
    exit 1
fi

REFERENCE=$1
PREFIX=$2

echo "Indexing $REFERENCE with Bowtie2..."
bowtie2-build "$REFERENCE" "$PREFIX"
echo "✅ Bowtie2 indexing complete!"
