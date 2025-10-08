#!/usr/bin/env python3
import scanpy as sc
import argparse, os, sys

parser = argparse.ArgumentParser(description="Run PCA on single or multiple h5ad files")
parser.add_argument("-i", "--input", nargs="+", required=True, help="Input .h5ad files")
parser.add_argument("-o", "--output_dir", required=True, help="Output directory")
args = parser.parse_args()

os.makedirs(args.output_dir, exist_ok=True)

for infile in args.input:
    try:
        adata = sc.read_h5ad(infile)
        sc.tl.pca(adata, svd_solver='arpack')
        outfile = os.path.join(args.output_dir, os.path.basename(infile))
        adata.write_h5ad(outfile)
        print(f"✅ PCA done for {infile} → {outfile}")
    except Exception as e:
        print(f"❌ Error processing {infile}: {e}", file=sys.stderr)
