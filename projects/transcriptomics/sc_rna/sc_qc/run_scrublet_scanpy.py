#!/usr/bin/env python3
"""
run_scrublet_scanpy.py
Professional doublet detection helper using Scrublet (Scanpy/AnnData compatible).

Creates per-sample outputs:
  - <sample>_scrublet.h5ad
  - <sample>_doublets.csv
"""
from __future__ import annotations
import argparse, os, sys, time, traceback

def log(msg, verbose=True):
    now = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{now}] {msg}", flush=True)

def check_packages():
    ok = True
    try: import scanpy as sc
    except Exception: log("ERROR: scanpy not installed. Install with `pip install scanpy`.", True); ok=False
    try: import scrublet as scr
    except Exception: log("ERROR: scrublet not installed. Install with `pip install scrublet`.", True); ok=False
    try: import numpy as np
    except Exception: log("ERROR: numpy not installed.", True); ok=False
    if not ok: sys.exit(1)

def safe_read_h5ad(path):
    import scanpy as sc
    if not os.path.exists(path): raise FileNotFoundError(path)
    try: return sc.read_h5ad(path)
    except Exception as e: raise RuntimeError(f"Failed to read {path}: {e}")

def ensure_counts_matrix_int(adata):
    import numpy as np, scipy.sparse as sp
    X = adata.X
    if sp.issparse(X): X = X.toarray()
    if not np.issubdtype(X.dtype, np.integer):
        log("Counts matrix not integer — rounding to nearest integer.")
        X = np.rint(X).astype(int)
    return X

def run_scrublet_for_file(fpath, outdir, rate=None, verbose=False):
    import numpy as np, pandas as pd, scanpy as sc, scrublet as scr
    sample = os.path.splitext(os.path.basename(fpath))[0]
    log(f"Processing {sample}", verbose)
    adata = safe_read_h5ad(fpath)
    X = ensure_counts_matrix_int(adata)
    scrub = scr.Scrublet(X, expected_doublet_rate=(rate if rate else 0.06))
    scores, preds = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)
    adata.obs["doublet_score"] = scores
    adata.obs["predicted_doublet"] = preds.astype(bool)
    os.makedirs(outdir, exist_ok=True)
    out_h5ad = os.path.join(outdir, f"{sample}_scrublet.h5ad")
    out_csv = os.path.join(outdir, f"{sample}_doublets.csv")
    adata.write(out_h5ad)
    adata.obs[["doublet_score","predicted_doublet"]].to_csv(out_csv)
    log(f"Wrote {out_h5ad} and {out_csv}", verbose)
    return out_h5ad, out_csv

def find_h5ad_files(d): return sorted([os.path.join(d,f) for f in os.listdir(d) if f.endswith(".h5ad")])

def parse_args():
    p = argparse.ArgumentParser(description="Run Scrublet on Scanpy h5ad files (single or directory mode).")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("-i","--input", help="Single input .h5ad")
    g.add_argument("-d","--dir", help="Directory containing .h5ad files")
    p.add_argument("-o","--outdir", required=True, help="Output directory")
    p.add_argument("-r","--rate", type=float, default=None, help="Expected doublet rate (e.g. 0.06)")
    p.add_argument("-v","--verbose", action="store_true", help="Verbose logging")
    return p.parse_args()

def main():
    args = parse_args(); check_packages(); os.makedirs(args.outdir, exist_ok=True)
    files = [args.input] if args.input else find_h5ad_files(args.dir)
    if not files: log("No .h5ad files found.", True); sys.exit(1)
    ok, fail = [], []
    for f in files:
        try: ok.append(run_scrublet_for_file(f, args.outdir, rate=args.rate, verbose=args.verbose))
        except Exception as e: fail.append((f,str(e))); log(f"ERROR on {f}: {e}"); traceback.print_exc()
    log(f"Done. Success: {len(ok)}; Failed: {len(fail)}")

if __name__ == "__main__": main()
