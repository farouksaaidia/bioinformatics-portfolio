#!/usr/bin/env bash
# Compress or decompress with gzip
# Usage: ./compress_decompress_gzip.sh [compress|decompress] file

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 compress|decompress file"
    exit 1
fi

if [ "$1" == "compress" ]; then
    gzip "$2"
    echo "✅ Compressed: $2.gz"
elif [ "$1" == "decompress" ]; then
    gunzip "$2"
    echo "✅ Decompressed: ${2%.gz}"
else
    echo "❌ Invalid option. Use compress|decompress."
    exit 1
fi
