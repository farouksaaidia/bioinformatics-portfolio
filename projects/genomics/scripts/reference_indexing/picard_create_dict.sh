#!/bin/bash
# Script: picard_create_dict.sh
# Description: Create a sequence dictionary (.dict) for GATK/Picard
# Usage: ./picard_create_dict.sh reference.fasta

if [ $# -ne 1 ]; then
    echo "Usage: $0 reference.fasta"
    exit 1
fi

REFERENCE=$1
OUTPUT=${REFERENCE%.fasta}.dict

echo "Creating dictionary for $REFERENCE with Picard..."
picard CreateSequenceDictionary R="$REFERENCE" O="$OUTPUT"
echo "✅ Picard dictionary created: $OUTPUT"
