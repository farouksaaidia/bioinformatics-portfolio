#!/usr/bin/env bash
# Index a bgzipped file with tabix
# Usage: ./index_tabix.sh input.vcf.gz

set -euo pipefail

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 input.bgz"
    exit 1
fi

tabix -p vcf "$1"
echo "✅ Tabix index created: $1.tbi"
