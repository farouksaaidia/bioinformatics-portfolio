#!/usr/bin/env python3
"""
Detect potential doublets in spatial transcriptomics data using Scrublet and spatial heuristics.
"""
import argparse, scanpy as sc, scrublet as scr, pandas as pd, numpy as np, os, sys

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input", required=True, help="Input .h5ad file")
parser.add_argument("-o","--output", required=True, help="Output .h5ad file")
parser.add_argument("--expected_doublet_rate", type=float, default=0.05, help="Expected doublet rate")
args = parser.parse_args()

if not os.path.exists(args.input):
    sys.exit(f"❌ Input file not found: {args.input}")

print("📥 Loading AnnData...")
adata = sc.read_h5ad(args.input)

# Run Scrublet
print("🧬 Running Scrublet doublet detection...")
counts_matrix = adata.X if isinstance(adata.X, np.ndarray) else adata.X.toarray()
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=args.expected_doublet_rate)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
adata.obs["doublet_score"] = doublet_scores
adata.obs["predicted_doublet"] = predicted_doublets

# Optional: spatial heuristic (density-based)
if "spatial" in adata.obsm.keys():
    coords = adata.obsm["spatial"]
    from sklearn.neighbors import NearestNeighbors
    nbrs = NearestNeighbors(n_neighbors=6).fit(coords)
    local_density = 1 / (np.mean(nbrs.kneighbors()[0], axis=1) + 1e-6)
    adata.obs["local_density"] = local_density
    adata.obs["density_adjusted_doublet_score"] = adata.obs["doublet_score"] * (local_density / np.median(local_density))

print(f"💾 Writing updated AnnData to {args.output}")
adata.write_h5ad(args.output)
print("✅ Doublet detection complete")
