#!/usr/bin/env python3
import scanpy as sc
import argparse, os, sys

parser = argparse.ArgumentParser(description="Run UMAP embedding on h5ad files with PCA done")
parser.add_argument("-i", "--input", nargs="+", required=True, help="Input .h5ad files with PCA")
parser.add_argument("-o", "--output_dir", required=True, help="Output directory for UMAP results")
args = parser.parse_args()

os.makedirs(args.output_dir, exist_ok=True)

for infile in args.input:
    try:
        adata = sc.read_h5ad(infile)
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
        sc.tl.umap(adata)
        outfile = os.path.join(args.output_dir, os.path.basename(infile))
        adata.write_h5ad(outfile)
        print(f"✅ UMAP done for {infile} → {outfile}")
    except Exception as e:
        print(f"❌ Error processing {infile}: {e}", file=sys.stderr)
