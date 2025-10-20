#!/usr/bin/env python3
"""
Perform batch correction/integration in Python using Scanpy + optional BBKNN or Scanorama.
Produces an AnnData with integrated PCA/UMAP.
"""
import scanpy as sc
import argparse, os, sys

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input", required=True, help="Input .h5ad (merged or single file with batch key)")
parser.add_argument("-o","--output", required=True, help="Output .h5ad (integrated)")
parser.add_argument("--batch_key", default="sample_id", help="obs column with batch labels (default: sample_id)")
parser.add_argument("--method", choices=["bbknn","scanorama","none"], default="bbknn", help="Integration method")
parser.add_argument("--n_pcs", type=int, default=30, help="Number of PCs")
args = parser.parse_args()

if not os.path.exists(args.input):
    sys.exit("❌ Input file not found")

adata = sc.read_h5ad(args.input)
if args.batch_key not in adata.obs:
    sys.exit(f"❌ Batch key '{args.batch_key}' not found in AnnData.obs")

print("ℹ️ Running PCA")
sc.pp.pca(adata, n_comps=args.n_pcs)

if args.method == "bbknn":
    try:
        import bbknn
    except Exception:
        sys.exit("❌ bbknn not installed. Install with `pip install bbknn`")
    print("🌀 Running BBKNN")
    bbknn.bbknn(adata, batch_key=args.batch_key)
    sc.tl.umap(adata)
elif args.method == "scanorama":
    try:
        import scanorama
    except Exception:
        sys.exit("❌ scanorama not installed. Install with `pip install scanorama`")
    print("🌀 Running Scanorama integration")
    sc.pp.neighbors(adata, n_pcs=args.n_pcs)
    sc.external.pp.scanorama_integrate(adata, key=args.batch_key)
    sc.tl.umap(adata)
else:
    print("⚠️ Skipping batch correction; computing neighbors on raw PCA")
    sc.pp.neighbors(adata, n_pcs=args.n_pcs)
    sc.tl.umap(adata)

adata.write_h5ad(args.output)
print("✅ Integration output written to", args.output)
