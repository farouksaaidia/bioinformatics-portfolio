#!/usr/bin/env python3
"""
Compare two annotation CSVs and compute metrics.
Assumes both CSVs have a shared cell identifier column (default: cell_id) and a label column (default: cell_type).
Outputs:
- {output}_confusion.csv
- {output}_metrics.txt (ARI, AMI, counts)
"""
import pandas as pd
import argparse
import sys
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score
import os

parser = argparse.ArgumentParser(description="Compare cell-type annotations between methods.")
parser.add_argument("-a", "--annot1", required=True, help="CSV with cluster annotations from method 1")
parser.add_argument("-b", "--annot2", required=True, help="CSV with cluster annotations from method 2")
parser.add_argument("-o", "--output", required=True, help="Output comparison results file prefix")
parser.add_argument("--id_col", default="cell_id", help="Column name for cell IDs (default: cell_id)")
parser.add_argument("--label_col", default="cell_type", help="Column name for labels (default: cell_type)")
args = parser.parse_args()

for f in [args.annot1, args.annot2]:
    if not os.path.exists(f):
        sys.exit(f"❌ Input file {f} not found")

try:
    ann1 = pd.read_csv(args.annot1)
    ann2 = pd.read_csv(args.annot2)
except Exception as e:
    sys.exit(f"❌ Failed to load input files: {e}")

if args.id_col not in ann1.columns or args.id_col not in ann2.columns:
    sys.exit(f"❌ id_col '{args.id_col}' not found in both files")

if args.label_col not in ann1.columns or args.label_col not in ann2.columns:
    sys.exit(f"❌ label_col '{args.label_col}' not found in both files")

merged = pd.merge(ann1[[args.id_col, args.label_col]], ann2[[args.id_col, args.label_col]],
                  on=args.id_col, suffixes=("_m1", "_m2"))
if merged.empty:
    sys.exit("❌ No matching cell IDs found between files after merge")

# compute metrics
ari = adjusted_rand_score(merged[f"{args.label_col}_m1"], merged[f"{args.label_col}_m2"])
ami = adjusted_mutual_info_score(merged[f"{args.label_col}_m1"], merged[f"{args.label_col}_m2"])

confmat = pd.crosstab(merged[f"{args.label_col}_m1"], merged[f"{args.label_col}_m2"])

confmat.to_csv(f"{args.output}_confusion.csv")

with open(f"{args.output}_metrics.txt", "w") as f:
    f.write(f"Adjusted Rand Index (ARI): {ari:.6f}\n")
    f.write(f"Adjusted Mutual Information (AMI): {ami:.6f}\n")
    f.write(f"Total cells compared: {len(merged)}\n")
    f.write("\n# Label counts method1\n")
    f.write(merged[f"{args.label_col}_m1"].value_counts().to_string())
    f.write("\n\n# Label counts method2\n")
    f.write(merged[f"{args.label_col}_m2"].value_counts().to_string())

print(f"✅ Comparison complete. ARI={ari:.6f}, AMI={ami:.6f}. Outputs: {args.output}_confusion.csv and {args.output}_metrics.txt")
