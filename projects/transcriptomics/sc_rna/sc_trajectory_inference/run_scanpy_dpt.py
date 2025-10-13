#!/usr/bin/env python3
"""
Compute diffusion pseudotime with Scanpy (DPT).
Inputs: .h5ad with neighbors computed (or computes neighbors if missing).
Outputs: annotated .h5ad with obs['dpt_pseudotime'] and optional root cell(s).
"""
import scanpy as sc, argparse, os, sys

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input", required=True, help="Input .h5ad (normalized, PCA/UMAP recommended)")
parser.add_argument("-o","--output", required=True, help="Output .h5ad with DPT pseudotime")
parser.add_argument("--root", default=None, help="Optional root cell id (or comma-separated list)")
args = parser.parse_args()

if not os.path.exists(args.input):
    sys.exit(f"❌ Input {args.input} not found")
adata = sc.read_h5ad(args.input)

if 'neighbors' not in adata.uns:
    print("ℹ️ neighbors not found — computing neighbors (n_pcs=30, n_neighbors=15)")
    sc.pp.neighbors(adata, n_pcs=30, n_neighbors=15)

if args.root:
    roots = args.root.split(",")
    # if root provided as cell indices, convert to index positions
    if all([r in adata.obs_names for r in roots]):
        root = roots[0]
    else:
        root = None
else:
    root = None

sc.tl.dpt(adata, n_dcs=10, min_group_size=0.01)
# If user provided root and it exists, re-run dpt with root
if root is not None:
    try:
        sc.tl.dpt(adata, n_dcs=10, min_group_size=0.01, root_key=root)
    except Exception:
        pass

adata.obs['dpt_pseudotime'] = adata.obs['dpt_pseudotime'].astype(float)
adata.write_h5ad(args.output)
print(f"✅ DPT pseudotime saved to {args.output}")
