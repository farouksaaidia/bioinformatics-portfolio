#!/usr/bin/env python3
import argparse, sys, pandas as pd, os

parser = argparse.ArgumentParser(description="Validate short-read RNA-seq sample sheet.")
parser.add_argument("-i","--input", required=True, help="Input sample sheet TSV (columns: sample_id, fastq1, [fastq2])")
args = parser.parse_args()

if not os.path.exists(args.input):
    sys.exit(f"❌ Sample sheet {args.input} not found.")

df = pd.read_csv(args.input, sep="\t")
required = ["sample_id", "fastq1"]
for col in required:
    if col not in df.columns:
        sys.exit(f"❌ Missing required column: {col}")

for idx,row in df.iterrows():
    if not os.path.exists(row.fastq1):
        sys.exit(f"❌ FASTQ not found: {row.fastq1}")
    if "fastq2" in df.columns and pd.notna(row.fastq2) and not os.path.exists(row.fastq2):
        sys.exit(f"❌ Paired FASTQ not found: {row.fastq2}")

print(f"✅ Valid sample sheet: {len(df)} samples validated.")
