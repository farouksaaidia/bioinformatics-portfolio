#!/usr/bin/env python3
"""
Automated cell-type annotation using CellTypist.
Supports model names (bundled) or local model files.

Outputs:
- annotated .h5ad with predictions in .obs['predicted_labels'] (or custom obs_key)
- optional CSV of predictions (use --pred_csv)
"""
import scanpy as sc
import celltypist
import argparse
import os
import sys
import pandas as pd

parser = argparse.ArgumentParser(description="Automated cell-type annotation using CellTypist.")
parser.add_argument("-i", "--input", required=True, help="Input clustered .h5ad file")
parser.add_argument("-m", "--model", default="Immune_All_Low.pkl", help="CellTypist model name or path (default: Immune_All_Low.pkl)")
parser.add_argument("-o", "--output", required=True, help="Output annotated .h5ad file")
parser.add_argument("--obs_key", default="predicted_labels", help="Observation key to store predictions (default: predicted_labels)")
parser.add_argument("--pred_csv", default=None, help="Optional path to save predictions CSV")
args = parser.parse_args()

if not os.path.exists(args.input):
    sys.exit(f"❌ Input file {args.input} not found")

print("📥 Loading data...")
adata = sc.read_h5ad(args.input)

# ensure adata has obs index
if adata.obs_names is None or len(adata.obs_names) == 0:
    sys.exit("❌ AnnData object has no obs_names (cell IDs).")

# determine model input
model_arg = args.model
if os.path.exists(model_arg):
    model_path = model_arg
    print(f"📦 Using local CellTypist model at {model_path}")
else:
    model_path = model_arg  # celltypist can accept known model names; celltypist will attempt to download if needed
    print(f"📦 Using CellTypist model name: {model_path} (if not available will attempt to download)")

print("🧠 Running CellTypist annotate... this may take a while for large datasets.")
try:
    predictions = celltypist.annotate(adata, model=model_path)
except Exception as e:
    sys.exit(f"❌ CellTypist annotation failed: {e}")

# Extract predicted labels robustly
pred_series = None
try:
    # celltypist.annotate may return a pandas.DataFrame with index matching cells
    if isinstance(predictions, pd.DataFrame):
        # Look for standard column names
        if "predicted_labels" in predictions.columns:
            pred_series = predictions["predicted_labels"]
        elif "labels" in predictions.columns:
            pred_series = predictions["labels"]
        else:
            # take first column as fallback
            pred_series = predictions.iloc[:, 0]
    else:
        # attempt attribute access
        pred_series = getattr(predictions, "predicted_labels", None)
        if pred_series is None:
            # some returns are Series-like
            pred_series = pd.Series(predictions, index=adata.obs_names)
except Exception as e:
    sys.exit(f"❌ Failed to extract predictions: {e}")

if pred_series is None:
    sys.exit("❌ Could not extract predicted labels from CellTypist result.")

# align to adata.obs_names
if not all(pred_series.index == adata.obs_names):
    # reindex if indexes differ
    pred_series = pred_series.reindex(adata.obs_names)
adata.obs[args.obs_key] = pred_series.values

# optional export CSV
if args.pred_csv:
    df = pd.DataFrame({ "cell_id": adata.obs_names, args.obs_key: adata.obs[args.obs_key].values })
    df.to_csv(args.pred_csv, index=False)
    print(f"💾 Predictions saved to {args.pred_csv}")

print(f"💾 Saving annotated object to {args.output}")
adata.write_h5ad(args.output)
print("✅ Automated annotation complete!")
