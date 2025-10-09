#!/usr/bin/env python3
"""
run_clustering_scanpy.py
Scanpy clustering wrapper (Leiden or Louvain). Supports single file or multiple inputs.

Usage examples:
  python run_clustering_scanpy.py -i sample.h5ad -o sample_clustered.h5ad -m leiden -r 0.5
  python run_clustering_scanpy.py -d input_dir/ -o out_dir/ -m louvain -r 0.4,0.6
"""
import argparse, os, sys
import scanpy as sc
import anndata as ad

def log(*args):
    print("[{}]".format(__import__("time").strftime("%Y-%m-%d %H:%M:%S")), *args, flush=True)

parser = argparse.ArgumentParser()
g = parser.add_mutually_exclusive_group(required=True)
g.add_argument("-i","--input", help="Input .h5ad file (single)")
g.add_argument("-d","--dir", help="Input directory with .h5ad files (batch)")
parser.add_argument("-o","--output", required=True, help="Output file (single) or output directory (batch)")
parser.add_argument("-m","--method", choices=["leiden","louvain"], default="leiden", help="Clustering method")
parser.add_argument("-r","--resolutions", default="0.5", help="Comma-separated resolution(s)")
parser.add_argument("-n","--npcs", type=int, default=30, help="Number of PCs to use")
args = parser.parse_args()

resolutions = [float(x) for x in args.resolutions.split(",")]

def process(infile, outpath):
    log("Reading", infile)
    if not os.path.exists(infile): raise FileNotFoundError(infile)
    adata = sc.read_h5ad(infile)
    if "X_pca" not in adata.obsm:
        log("PCA missing — running PCA with", args.npcs, "components")
        sc.tl.pca(adata, n_comps=args.npcs)
    sc.pp.neighbors(adata, n_pcs=args.npcs)
    for r in resolutions:
        if args.method == "leiden":
            sc.tl.leiden(adata, resolution=r, key_added=f"leiden_r{r}")
            adata.obs[f"cluster_resolution_{r}"] = adata.obs[f"leiden_r{r}"]
        else:
            sc.tl.louvain(adata, resolution=r, key_added=f"louvain_r{r}")
            adata.obs[f"cluster_resolution_{r}"] = adata.obs[f"louvain_r{r}"]
    # set final cluster column to last resolution
    final_col = f"cluster_resolution_{resolutions[-1]}"
    adata.obs["seurat_clusters"] = adata.obs[final_col].astype(str)
    os.makedirs(os.path.dirname(outpath) if os.path.dirname(outpath) else ".", exist_ok=True)
    adata.write_h5ad(outpath)
    log("Wrote clustered file to", outpath)

try:
    if args.input:
        outfile = args.output if args.output.endswith(".h5ad") else args.output
        process(args.input, outfile)
    else:
        if not os.path.isdir(args.dir): raise NotADirectoryError(args.dir)
        os.makedirs(args.output, exist_ok=True)
        files = [os.path.join(args.dir,f) for f in os.listdir(args.dir) if f.endswith(".h5ad")]
        if not files: raise FileNotFoundError("No .h5ad files in input directory")
        for f in files:
            outfile = os.path.join(args.output, os.path.basename(f))
            process(f, outfile)
except Exception as e:
    log("ERROR:", e)
    sys.exit(1)
