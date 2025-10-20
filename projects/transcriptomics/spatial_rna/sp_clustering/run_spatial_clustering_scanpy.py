#!/usr/bin/env python3
"""
Perform spatial clustering on AnnData object using Scanpy (Leiden/Louvain).
"""
import scanpy as sc
import argparse, os, sys

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input", required=True, help="Input .h5ad file")
parser.add_argument("-o","--output", required=True, help="Output clustered .h5ad file")
parser.add_argument("--method", choices=["leiden","louvain"], default="leiden", help="Clustering method")
parser.add_argument("--resolution", type=float, default=0.6, help="Clustering resolution")
args = parser.parse_args()

if not os.path.exists(args.input):
    sys.exit(f"❌ Input file not found: {args.input}")

print("📥 Loading data...")
adata = sc.read_h5ad(args.input)

if 'neighbors' not in adata.uns:
    print("ℹ️ No neighbors found — computing with default params.")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)

print(f"🧭 Running {args.method} clustering at resolution {args.resolution}")
if args.method == "leiden":
    sc.tl.leiden(adata, resolution=args.resolution, key_added="spatial_clusters")
else:
    sc.tl.louvain(adata, resolution=args.resolution, key_added="spatial_clusters")

print("💾 Saving clustered AnnData")
adata.write_h5ad(args.output)
print("✅ Clustering complete.")
