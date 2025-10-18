#!/usr/bin/env python3
"""
Load Visium / 10X spatial outputs into an AnnData .h5ad.
- Supports pointing at sample folder or sample/outs folder.
- Preserves spatial_coords and image metadata when available.
"""
import scanpy as sc
import anndata
import argparse
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-s","--sample_dir", required=True, help="Visium sample directory or 'outs' path")
parser.add_argument("-n","--sample_name", required=True, help="Sample name")
parser.add_argument("-o","--output", required=True, help="Output .h5ad path")
args = parser.parse_args()

sample_dir = args.sample_dir
if not os.path.exists(sample_dir):
    # try sample_dir/outs
    alt = os.path.join(sample_dir, "outs")
    if os.path.exists(alt):
        sample_dir = alt
    else:
        sys.exit(f"❌ Sample path not found: {args.sample_dir}")

print(f"📥 Loading Visium data from {sample_dir}")
try:
    adata = sc.read_visium(sample_dir)
except Exception as e:
    sys.exit(f"❌ Scanpy read_visium failed: {e}")

# Add sample name
adata.obs['sample_id'] = args.sample_name

# Basic QC: percent mt (if MT- genes present)
mt_genes = [g for g in adata.var_names if g.upper().startswith("MT-") or g.upper().startswith("MT")]
if len(mt_genes) > 0:
    adata.obs['percent_mt'] = (adata[:, mt_genes].X.sum(axis=1).A1 if hasattr(adata.X, 'A1') else adata[:, mt_genes].X.sum(axis=1)) / adata.X.sum(axis=1)

print(f"💾 Writing AnnData to {args.output}")
adata.write_h5ad(args.output)
print("✅ Done.")
