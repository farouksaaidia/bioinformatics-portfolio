#!/usr/bin/env python3
"""
Compute PCA, neighbors, and UMAP for spatial AnnData objects (Scanpy).
Preserves spatial coordinates in .obsm['spatial'].
"""
import scanpy as sc
import argparse
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input", required=True, help="Input .h5ad (normalized recommended)")
parser.add_argument("-o","--output", required=True, help="Output .h5ad with UMAP")
parser.add_argument("--n_pcs", type=int, default=30, help="Number of PCs to compute")
parser.add_argument("--n_neighbors", type=int, default=15, help="n_neighbors for neighbor graph")
parser.add_argument("--min_dist", type=float, default=0.3, help="min_dist for UMAP")
args = parser.parse_args()

if not os.path.exists(args.input):
    sys.exit(f"❌ Input not found: {args.input}")

print("📥 Loading AnnData...")
adata = sc.read_h5ad(args.input)

# PCA
if 'X_pca' not in adata.obsm_keys():
    print(f"ℹ️ Running PCA (n_comps={args.n_pcs})")
    sc.pp.pca(adata, n_comps=args.n_pcs, svd_solver='arpack')
else:
    print("ℹ️ PCA already present in obsm['X_pca'] — skipping PCA")

print(f"ℹ️ Computing neighbors (n_neighbors={args.n_neighbors})")
sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, n_pcs=args.n_pcs)

print("🧭 Computing UMAP")
sc.tl.umap(adata, min_dist=args.min_dist)

print(f"💾 Writing AnnData to {args.output}")
adata.write_h5ad(args.output)
print("✅ UMAP saved.")
