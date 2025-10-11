#!/usr/bin/env python3
"""
Automated cell-type annotation using CellTypist (and fallback-friendly handling).
Requirements: scanpy, celltypist, pandas
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
parser.add_argument("--pred-col", default="predicted_labels", help="Column name to store predictions in adata.obs")
parser.add_argument("--prob-col", default="prediction_prob", help="Optional column name for prediction confidence (if available)")
args = parser.parse_args()

if not os.path.exists(args.input):
    sys.exit(f"❌ Input file {args.input} not found")

# load AnnData
print("📥 Loading data...")
try:
    adata = sc.read_h5ad(args.input)
except Exception as e:
    sys.exit(f"❌ Failed to read input .h5ad: {e}")

# ensure obs exists
if adata.n_obs == 0:
    sys.exit("❌ Input AnnData contains no cells")

# Try to load model: if args.model is a path, check exists; otherwise use as model name for celltypist
model_arg = args.model
if os.path.exists(model_arg):
    model_path = model_arg
    print(f"🗂 Using local model at {model_path}")
else:
    # try to fetch/confirm model name via celltypist.models
    model_path = model_arg
    print(f"🧠 Using CellTypist model identifier: {model_arg} (celltypist will attempt to resolve/download)")

print("🔎 Running CellTypist annotation (this may take a moment)...")
try:
    # celltypist.annotate returns a pandas DataFrame-like or object depending on versions.
    predictions = celltypist.annotate(adata, model=model_path, majority_voting=True, return_prob=True)
except Exception as e:
    # try a more basic call without return_prob
    try:
        predictions = celltypist.annotate(adata, model=model_path, majority_voting=True)
    except Exception as e2:
        sys.exit(f"❌ CellTypist annotation failed: {e2}")

# Extract predicted labels robustly
pred_series = None
prob_series = None

# If predictions is a DataFrame (newer versions), expect 'predicted_labels' or 'label' column
if isinstance(predictions, pd.DataFrame):
    if 'predicted_labels' in predictions.columns:
        pred_series = predictions['predicted_labels'].astype(str)
    elif 'label' in predictions.columns:
        pred_series = predictions['label'].astype(str)
    elif predictions.shape[1] >= 1:
        pred_series = predictions.iloc[:, 0].astype(str)
    # try to find a probability/confidence column
    for c in ('probability', 'score', 'confidence', 'prob'):
        if c in predictions.columns:
            prob_series = predictions[c]
            break
# If predictions has attributes (older API)
else:
    # attempt to access attributes conservatively
    pred_series = getattr(predictions, 'predicted_labels', None)
    if pred_series is None:
        pred_series = getattr(predictions, 'labels', None)
    if pred_series is None and hasattr(predictions, '__iter__'):
        try:
            # try converting to Series
            pred_series = pd.Series(list(predictions)).astype(str)
        except Exception:
            pred_series = None
    # confidence/prob
    prob_series = getattr(predictions, 'probability', None) or getattr(predictions, 'probabilities', None)

if pred_series is None:
    sys.exit("❌ Could not extract predicted labels from CellTypist output (incompatible API).")

# align predictions to adata.obs index if possible
if len(pred_series) == adata.n_obs:
    adata.obs[args.pred_col] = pred_series.values
else:
    # try if predictions index matches adata.obs_names
    try:
        pred_series.index = pred_series.index.astype(str)
        obs_index = pd.Index(map(str, adata.obs_names))
        aligned = pred_series.reindex(obs_index)
        if aligned.isnull().all():
            # fallback: store as-is with warning
            adata.obs[args.pred_col] = pd.Series(pred_series.values[:adata.n_obs], index=adata.obs_names)
            print("⚠️ Warning: prediction length did not match number of cells; truncated/filled to fit AnnData.")
        else:
            adata.obs[args.pred_col] = aligned.values
    except Exception:
        adata.obs[args.pred_col] = pd.Series(pred_series.values[:adata.n_obs], index=adata.obs_names)
        print("⚠️ Warning: prediction length did not match number of cells; truncated/filled to fit AnnData.")

if prob_series is not None:
    # attempt similar alignment for probability/confidence
    if len(prob_series) == adata.n_obs:
        adata.obs[args.prob_col] = prob_series.values
    else:
        try:
            alignedp = prob_series.reindex(pd.Index(map(str, adata.obs_names)))
            adata.obs[args.prob_col] = alignedp.values
        except Exception:
            print("⚠️ Could not align probability/confidence vector to AnnData; skipping.")

print(f"💾 Saving annotated object to {args.output}")
try:
    adata.write_h5ad(args.output)
except Exception as e:
    sys.exit(f"❌ Failed to write output .h5ad: {e}")

print("✅ Automated annotation complete!")
