#!/usr/bin/env python3
"""
Compare cell-type annotations between two methods.
Outputs:
 - <output>_confusion.csv (cross-tabulation)
 - <output>_metrics.txt (ARI and basic stats)
Requirements: pandas, scikit-learn
"""
import pandas as pd
import argparse
import sys
import os
from sklearn.metrics import adjusted_rand_score

parser = argparse.ArgumentParser(description="Compare cell-type annotations between methods.")
parser.add_argument("-a", "--annot1", required=True, help="CSV with annotations from method 1 (must include cell_id and cell_type columns unless --id-col/--type-col provided)")
parser.add_argument("-b", "--annot2", required=True, help="CSV with annotations from method 2")
parser.add_argument("-o", "--output", required=True, help="Output comparison results file prefix")
parser.add_argument("--id-col", default="cell_id", help="Column name for cell identifier (default: cell_id)")
parser.add_argument("--type-col", default="cell_type", help="Column name for annotation label (default: cell_type)")
args = parser.parse_args()

def load_csv(path):
    if not os.path.exists(path):
        sys.exit(f"❌ File not found: {path}")
    try:
        return pd.read_csv(path)
    except Exception as e:
        sys.exit(f"❌ Failed to read {path}: {e}")

ann1 = load_csv(args.annot1)
ann2 = load_csv(args.annot2)

if args.id_col not in ann1.columns or args.type_col not in ann1.columns:
    sys.exit(f"❌ annot1 must contain columns: {args.id_col}, {args.type_col}")
if args.id_col not in ann2.columns or args.type_col not in ann2.columns:
    sys.exit(f"❌ annot2 must contain columns: {args.id_col}, {args.type_col}")

# Merge on cell id
merged = pd.merge(ann1[[args.id_col, args.type_col]], ann2[[args.id_col, args.type_col]],
                  on=args.id_col, suffixes=("_m1", "_m2"), how="inner")

if merged.empty:
    sys.exit("❌ No overlapping cell IDs found between the two annotation files.")

# Compute ARI
ari = adjusted_rand_score(merged[f"{args.type_col}_m1"], merged[f"{args.type_col}_m2"])

# Confusion matrix (cross-tab)
confmat = pd.crosstab(merged[f"{args.type_col}_m1"], merged[f"{args.type_col}_m2"])

# Save outputs
conf_csv = f"{args.output}_confusion.csv"
metrics_txt = f"{args.output}_metrics.txt"

confmat.to_csv(conf_csv)
with open(metrics_txt, "w") as f:
    f.write(f"Adjusted Rand Index: {ari:.4f}\n")
    f.write(f"Method1 annotations: {merged[f'{args.type_col}_m1'].nunique()} classes\n")
    f.write(f"Method2 annotations: {merged[f'{args.type_col}_m2'].nunique()} classes\n")
    f.write(f"Compared cells (overlap): {len(merged)}\n")

print(f"✅ Comparison complete. ARI={ari:.4f}")
print(f"Outputs written: {conf_csv}, {metrics_txt}")
