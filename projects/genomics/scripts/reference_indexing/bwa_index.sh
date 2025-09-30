#!/bin/bash
# Script: bwa_index.sh
# Description: Index a reference FASTA for BWA alignment
# Usage: ./bwa_index.sh reference.fasta

if [ $# -ne 1 ]; then
    echo "Usage: $0 reference.fasta"
    exit 1
fi

REFERENCE=$1

echo "Indexing $REFERENCE with BWA..."
bwa index $REFERENCE
echo "BWA indexing complete!"
