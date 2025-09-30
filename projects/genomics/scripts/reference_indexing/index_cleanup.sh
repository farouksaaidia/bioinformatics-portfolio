#!/bin/bash
# Script: index_cleanup.sh
# Description: Remove old index files to avoid conflicts before re-indexing
# Usage: ./index_cleanup.sh reference.fasta

if [ $# -ne 1 ]; then
    echo "Usage: $0 reference.fasta"
    exit 1
fi

REFERENCE=$1
BASENAME=${REFERENCE%.fasta}

echo "🧹 Cleaning up index files for $REFERENCE..."
rm -f "$REFERENCE".fai "$BASENAME".dict "$BASENAME".*bwt "$BASENAME".*sa "$BASENAME".*bt2

echo "✅ Old index files removed!"
