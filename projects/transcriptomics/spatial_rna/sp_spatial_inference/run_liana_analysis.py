#!/usr/bin/env python3
"""
Run ligand–receptor inference using LIANA (multi-method consensus analysis).
"""
import liana as li
import scanpy as sc
import argparse, os, sys

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input", required=True, help="Input .h5ad file (with cell_type column)")
parser.add_argument("-o","--output", required=True, help="Output CSV file for ligand–receptor results")
args = parser.parse_args()

if not os.path.exists(args.input):
    sys.exit(f"❌ Input file not found: {args.input}")

print("📥 Loading AnnData...")
adata = sc.read_h5ad(args.input)

if "cell_type" not in adata.obs:
    sys.exit("❌ 'cell_type' column required in AnnData.obs")

print("🧩 Running LIANA consensus ligand–receptor inference...")
li.mt.rank_aggregate(adata, resource="consensus")

print("💾 Saving ligand–receptor results...")
adata.uns['liana_res'].to_csv(args.output, index=False)
print("✅ LIANA analysis complete ->", args.output)
