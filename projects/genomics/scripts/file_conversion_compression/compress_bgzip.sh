#!/usr/bin/env bash
# Compress a file using bgzip (block compression)
# Usage: ./compress_bgzip.sh input.vcf output.vcf.gz

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input_file output_file.gz"
    exit 1
fi

bgzip -c "$1" > "$2"
echo "✅ File compressed with bgzip: $2"
