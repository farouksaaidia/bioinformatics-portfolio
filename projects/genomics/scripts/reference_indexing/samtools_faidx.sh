#!/bin/bash
# Script: samtools_faidx.sh
# Description: Create FASTA index (.fai) using samtools
# Usage: ./samtools_faidx.sh reference.fasta

if [ $# -ne 1 ]; then
    echo "Usage: $0 reference.fasta"
    exit 1
fi

REFERENCE=$1

echo "Indexing $REFERENCE with samtools..."
samtools faidx "$REFERENCE"
echo "✅ Samtools FASTA indexing complete!"
