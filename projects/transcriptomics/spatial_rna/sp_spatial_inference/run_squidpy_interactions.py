#!/usr/bin/env python3
"""
Compute spatial co-occurrence and cell–cell proximity enrichment using Squidpy.
"""
import squidpy as sq
import scanpy as sc
import argparse, os, sys

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input", required=True, help="Input .h5ad file with spatial coords")
parser.add_argument("-o","--output_prefix", required=True, help="Output prefix for results")
args = parser.parse_args()

adata = sc.read_h5ad(args.input)
print("📥 Computing spatial neighbors...")
sq.gr.spatial_neighbors(adata)

print("🧩 Running spatial co-occurrence...")
sq.gr.co_occurrence(adata, cluster_key="cell_type")

print("💾 Saving co-occurrence results...")
adata.uns["co_occurrence"]["occ"].to_csv(f"{args.output_prefix}_cooccurrence.csv")
print("✅ Squidpy interaction analysis complete ->", args.output_prefix)
