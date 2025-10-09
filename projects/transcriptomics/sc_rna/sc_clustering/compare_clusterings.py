#!/usr/bin/env python3
"""
compare_clusterings.py
Compare two clustering columns in an AnnData or Seurat-exported CSV.
Computes ARI and outputs contingency table (CSV + optional heatmap).
Usage:
  python compare_clusterings.py -i input.h5ad -c1 leiden_r0.4 -c2 louvain_r0.6 -o out_prefix
"""
import argparse, os, sys
import scanpy as sc
import pandas as pd
from sklearn.metrics import adjusted_rand_score
import seaborn as sns
import matplotlib.pyplot as plt

def log(*args): print("[{}]".format(__import__("time").strftime("%Y-%m-%d %H:%M:%S")), *args)

p = argparse.ArgumentParser()
p.add_argument("-i","--input", required=True, help="Input .h5ad")
p.add_argument("-c1","--col1", required=True, help="Clustering column 1 in .obs")
p.add_argument("-c2","--col2", required=True, help="Clustering column 2 in .obs")
p.add_argument("-o","--out", required=True, help="Output prefix (folder or file prefix)")
args = p.parse_args()

if not os.path.exists(args.input): log("Input not found:", args.input); sys.exit(1)
adata = sc.read_h5ad(args.input)
if args.col1 not in adata.obs.columns or args.col2 not in adata.obs.columns:
    log("Clustering columns not found in adata.obs. Available columns:", list(adata.obs.columns))
    sys.exit(1)

labels1 = adata.obs[args.col1].astype(str)
labels2 = adata.obs[args.col2].astype(str)
ari = adjusted_rand_score(labels1, labels2)
log("Adjusted Rand Index (ARI):", ari)

ct = pd.crosstab(labels1, labels2)
out_folder = os.path.dirname(args.out) if os.path.dirname(args.out) else "."
os.makedirs(out_folder, exist_ok=True)
ct.to_csv(args.out + "_contingency.csv")
with plt.style.context('seaborn'):
    plt.figure(figsize=(8,6))
    sns.heatmap(ct, cmap="viridis", annot=False)
    plt.title(f"Contingency: {args.col1} vs {args.col2} (ARI={ari:.3f})")
    plt.tight_layout()
    plt.savefig(args.out + "_contingency_heatmap.png", dpi=150)
log("Wrote contingency CSV and heatmap with prefix:", args.out)
