#!/usr/bin/env python3
"""
generate_qc_report.py — Professional QC report generator for single-cell AnnData (.h5ad)

Features:
 - Single-file mode (-i) or directory mode (-d)
 - Generates per-sample PNG plots (violin: n_genes/n_counts/pct_mt; scatter; doublet hist)
 - Produces per-sample CSV summary (basic metrics)
 - Produces a combined HTML report (qc_report.html) embedding images and tables
 - Robust argument parsing, dependency checks, logging, error handling

Usage examples:
  Single: ./generate_qc_report.py -i sample_scrublet.h5ad -o ./qc_reports
  Batch : ./generate_qc_report.py -d ./h5ad_dir -o ./qc_reports --workers 4
"""
from __future__ import annotations
import argparse
import os
import sys
import time
import traceback
from typing import List, Dict

def log(msg: str, level: str = "INFO"):
    now = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{now}] {level}: {msg}", flush=True)

def check_dependencies():
    missing = []
    try:
        import scanpy as sc  # noqa: F401
    except Exception:
        missing.append("scanpy")
    try:
        import pandas as pd  # noqa: F401
    except Exception:
        missing.append("pandas")
    try:
        import matplotlib  # noqa: F401
    except Exception:
        missing.append("matplotlib")
    try:
        import seaborn  # noqa: F401
    except Exception:
        missing.append("seaborn")
    if missing:
        log("Missing Python packages: " + ", ".join(missing), "ERROR")
        log("Install them via: pip install scanpy pandas matplotlib seaborn", "ERROR")
        sys.exit(1)

def safe_read_h5ad(path: str):
    import scanpy as sc
    if not os.path.exists(path):
        raise FileNotFoundError(f"{path} not found")
    try:
        adata = sc.read_h5ad(path)
        return adata
    except Exception as e:
        raise RuntimeError(f"Failed to read {path} as h5ad: {e}")

def normalize_obs_names(adata):
    # create/normalize common metric names used later
    obs = adata.obs
    # detect n_genes-like columns
    for cand in ["n_genes", "n_genes_by_counts", "nFeature_RNA", "n_genes_by_counts"]:
        if cand in obs.columns:
            obs["n_genes_common"] = obs[cand]
            break
    if "n_genes_common" not in obs.columns:
        if hasattr(adata, 'n_vars'):
            obs["n_genes_common"] = adata.var.shape[0]  # fallback (not ideal)
    # counts
    for cand in ["total_counts", "n_counts", "nCount_RNA", "nCount_RNA"]:
        if cand in obs.columns:
            obs["n_counts_common"] = obs[cand]
            break
    # percent mt
    for cand in ["pct_counts_mt", "percent.mt", "pct_mt", "percent_mito"]:
        if cand in obs.columns:
            obs["pct_mt_common"] = obs[cand]
            break

def plot_violin(df, outpath, title_suffix=""):
    import matplotlib.pyplot as plt
    import seaborn as sns
    plt.figure(figsize=(10,4))
    cols = []
    if "n_genes_common" in df.columns: cols.append("n_genes_common")
    if "n_counts_common" in df.columns: cols.append("n_counts_common")
    if "pct_mt_common" in df.columns: cols.append("pct_mt_common")
    n = len(cols)
    if n == 0:
        log("No common QC columns found for violin plot", "WARN")
        return None
    fig, axes = plt.subplots(1, n, figsize=(5*n, 4))
    if n == 1:
        axes = [axes]
    for ax, col in zip(axes, cols):
        sns.violinplot(y=df[col], ax=ax, inner="quartile")
        ax.set_title(col.replace("_common","") + (" " + title_suffix if title_suffix else ""))
        if col == "n_counts_common":
            ax.set_yscale("log")
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    return outpath

def plot_scatter(df, outpath):
    import matplotlib.pyplot as plt
    import seaborn as sns
    if ("n_counts_common" not in df.columns) or ("n_genes_common" not in df.columns):
        log("Skipping scatter plot: required columns not found", "WARN")
        return None
    plt.figure(figsize=(6,5))
    sns.scatterplot(x=df["n_counts_common"], y=df["n_genes_common"], hue=df.get("pct_mt_common", None), palette="viridis", s=8)
    plt.xscale("log")
    plt.xlabel("n_counts")
    plt.ylabel("n_genes")
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()
    return outpath

def plot_doublet_hist(df, outpath):
    import matplotlib.pyplot as plt
    import seaborn as sns
    if "doublet_score" not in df.columns:
        log("No doublet_score in data; skipping histogram", "INFO")
        return None
    plt.figure(figsize=(6,3))
    sns.histplot(df["doublet_score"].dropna(), bins=50)
    plt.xlabel("doublet_score")
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()
    return outpath

def summarize_metrics(df):
    import numpy as np
    import pandas as pd
    d = {
        "n_cells": int(df.shape[0]),
        "median_n_genes": int(df["n_genes_common"].median()) if "n_genes_common" in df.columns else None,
        "median_n_counts": int(df["n_counts_common"].median()) if "n_counts_common" in df.columns else None,
        "median_pct_mt": float(df["pct_mt_common"].median()) if "pct_mt_common" in df.columns else None
    }
    return pd.DataFrame([d])

def make_html_report(sample_info: Dict[str, Dict], outdir: str, out_html="qc_report.html"):
    html_lines = []
    html_lines.append("<!doctype html>")
    html_lines.append("<html><head><meta charset='utf-8'><title>scRNA-seq QC Report</title>")
    html_lines.append("<style> body{font-family:Arial,sans-serif;margin:20px;} img{border:1px solid #ddd;padding:4px;margin:8px;} .sample{margin-bottom:40px;} table{border-collapse:collapse;} th,td{border:1px solid #ccc;padding:6px;} </style>")
    html_lines.append("</head><body>")
    html_lines.append(f"<h1>scRNA-seq QC Report</h1>")
    html_lines.append(f"<p>Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}</p>")
    for sample, info in sample_info.items():
        html_lines.append(f"<div class='sample'><h2>{sample}</h2>")
        if info.get("pngs"):
            for p in info["pngs"]:
                html_lines.append(f"<img src='{os.path.basename(p)}' style='width:420px;'/>")
        if info.get("summary_table") is not None:
            # embed basic table
            html_lines.append("<h3>Summary</h3>")
            html_lines.append(info["summary_table"].to_html(index=False))
        if info.get("csv"):
            html_lines.append(f"<p>CSV summary: <a href='{os.path.basename(info['csv'])}'>{os.path.basename(info['csv'])}</a></p>")
        html_lines.append("</div><hr/>")
    html_lines.append("</body></html>")
    outpath = os.path.join(outdir, out_html)
    with open(outpath, "w") as fh:
        fh.write("\n".join(html_lines))
    log(f"Wrote HTML report to {outpath}")
    return outpath

def process_file(fpath: str, outdir: str):
    import scanpy as sc
    import pandas as pd
    sample = os.path.splitext(os.path.basename(fpath))[0]
    log(f"Processing {sample}")
    adata = safe_read_h5ad(fpath)
    normalize_obs_names(adata)
    obs = adata.obs.copy()
    # ensure numeric columns exist
    # prefer existing common columns or recalc basic metrics
    if "n_genes_common" not in obs.columns:
        if hasattr(adata, 'X'):
            log("n_genes_common not found — computing per-cell gene counts from matrix.")
            obs["n_genes_common"] = (adata.X > 0).sum(axis=1).A1 if hasattr(adata.X, "A1") else (adata.X > 0).sum(axis=1)
    if "n_counts_common" not in obs.columns:
        if hasattr(adata, 'X'):
            log("n_counts_common not found — computing per-cell count sums from matrix.")
            obs["n_counts_common"] = adata.X.sum(axis=1).A1 if hasattr(adata.X, "A1") else adata.X.sum(axis=1)
    # standardize names
    if "pct_mt_common" not in obs.columns:
        # try derive if MT genes exist in var names
        mt_genes = [g for g in adata.var_names if g.upper().startswith("MT-") or g.upper().startswith("MT.")]
        if mt_genes:
            log("Computing percent mitochondrial (pct_mt_common) from MT gene set.")
            sc.pp.calculate_qc_metrics(adata, qc_vars=[mt_genes], inplace=True)
            # many versions put result in adata.obs['pct_counts_in_top_...'] — try common names
            for cand in adata.obs.columns:
                if "pct_counts" in cand and "mt" in cand.lower():
                    obs["pct_mt_common"] = adata.obs[cand]
                    break
    # convert to DataFrame numeric types where possible
    obs = obs.copy()
    # create output PNGs
    pngs = []
    violin_png = os.path.join(outdir, f"{sample}_violin.png")
    s = plot_violin(obs, violin_png, title_suffix=sample)
    if s: pngs.append(s)
    scatter_png = os.path.join(outdir, f"{sample}_scatter.png")
    s2 = plot_scatter(obs, scatter_png)
    if s2: pngs.append(s2)
    doublet_png = os.path.join(outdir, f"{sample}_doublet_hist.png")
    s3 = plot_doublet_hist(obs, doublet_png)
    if s3: pngs.append(s3)
    # summary CSV
    summary_df = summarize_metrics(obs)
    csv_path = os.path.join(outdir, f"{sample}_qc_summary.csv")
    summary_df.to_csv(csv_path, index=False)
    log(f"Wrote summary CSV to {csv_path}")
    return {"pngs": pngs, "summary_table": summary_df, "csv": csv_path}

def find_h5ad_files(dirpath: str) -> List[str]:
    return sorted([os.path.join(dirpath, f) for f in os.listdir(dirpath) if f.endswith(".h5ad")])

def parse_args():
    p = argparse.ArgumentParser(description="Generate QC report for scRNA-seq h5ad files (single or directory).")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("-i", "--input", help="Single input .h5ad file")
    g.add_argument("-d", "--dir", help="Directory containing .h5ad files")
    p.add_argument("-o", "--outdir", required=True, help="Output directory")
    p.add_argument("-w", "--workers", type=int, default=1, help="Number of parallel workers (currently serial, reserved)")
    return p.parse_args()

def main():
    check_dependencies()
    args = parse_args()
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
    files = []
    if args.input:
        if not os.path.isfile(args.input):
            log(f"Input file not found: {args.input}", "ERROR"); sys.exit(1)
        files = [args.input]
    else:
        if not os.path.isdir(args.dir):
            log(f"Input directory not found: {args.dir}", "ERROR"); sys.exit(1)
        files = find_h5ad_files(args.dir)
        if not files:
            log(f"No .h5ad files found in {args.dir}", "ERROR"); sys.exit(1)
    sample_info = {}
    for f in files:
        try:
            info = process_file(f, outdir)
            sample_name = os.path.splitext(os.path.basename(f))[0]
            sample_info[sample_name] = info
        except Exception as e:
            log(f"ERROR processing {f}: {e}", "ERROR")
            traceback.print_exc()
    # write combined HTML
    make_html_report(sample_info, outdir)
    log("All done.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
