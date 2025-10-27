#!/usr/bin/env bash
set -euo pipefail

usage() { echo "Usage: $0 -i <fast5_dir> -o <output.csv>"; exit 1; }

while getopts "i:o:" opt; do
  case ${opt} in
    i) INPUT=$OPTARG ;;
    o) OUTPUT=$OPTARG ;;
    *) usage ;;
  esac
done
[[ -z "$INPUT" || -z "$OUTPUT" ]] && usage

mkdir -p "$(dirname "$OUTPUT")"

echo "🧬 Extracting Poly(A) tail lengths..."

python3 - <<PY
import os, csv
from ont_fast5_api.fast5_interface import get_fast5_file

rows=[]
for root, _, files in os.walk("$INPUT"):
    for f in files:
        if not f.endswith(".fast5"): continue
        with get_fast5_file(os.path.join(root,f), mode="r") as h:
            try:
                ds = h.get_analysis_dataset("Basecall_1D_000", "BaseCalled_template")
                rows.append([f, ds.attrs.get("poa_length", None)])
            except Exception:
                pass
with open("$OUTPUT","w",newline="") as out:
    w=csv.writer(out); w.writerow(["read","polya_length"]); w.writerows(rows)
PY

echo "✅ Poly(A) summary → $OUTPUT"
