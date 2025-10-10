#!/usr/bin/env python3
import scanpy as sc
import celltypist
import argparse
import os
import sys

parser = argparse.ArgumentParser(description="Automated cell-type annotation using CellTypist/Azimuth reference models.")
parser.add_argument("-i", "--input", required=True, help="Input clustered .h5ad file")
parser.add_argument("-m", "--model", default="Immune_All_Low.pkl", help="CellTypist model name (default: Immune_All_Low.pkl)")
parser.add_argument("-o", "--output", required=True, help="Output annotated .h5ad file")
args = parser.parse_args()

if not os.path.exists(args.input):
    sys.exit(f"❌ Input file {args.input} not found")

print("📥 Loading data...")
adata = sc.read_h5ad(args.input)

print(f"🧠 Running CellTypist with model {args.model}...")
predictions = celltypist.annotate(adata, model=args.model)
adata.obs['predicted_labels'] = predictions.predicted_labels

print(f"💾 Saving annotated object to {args.output}")
adata.write_h5ad(args.output)
print("✅ Automated annotation complete!")
