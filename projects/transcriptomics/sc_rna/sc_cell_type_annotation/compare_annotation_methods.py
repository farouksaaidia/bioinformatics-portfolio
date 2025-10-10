#!/usr/bin/env python3
import pandas as pd
import argparse
import sys
from sklearn.metrics import adjusted_rand_score, confusion_matrix

parser = argparse.ArgumentParser(description="Compare cell-type annotations between methods.")
parser.add_argument("-a", "--annot1", required=True, help="CSV with cluster annotations from method 1")
parser.add_argument("-b", "--annot2", required=True, help="CSV with cluster annotations from method 2")
parser.add_argument("-o", "--output", required=True, help="Output comparison results file prefix")
args = parser.parse_args()

try:
    ann1 = pd.read_csv(args.annot1)
    ann2 = pd.read_csv(args.annot2)
except Exception as e:
    sys.exit(f"❌ Failed to load input files: {e}")

merged = pd.merge(ann1, ann2, on="cell_id", suffixes=("_m1", "_m2"))
ari = adjusted_rand_score(merged.cell_type_m1, merged.cell_type_m2)
confmat = pd.crosstab(merged.cell_type_m1, merged.cell_type_m2)

confmat.to_csv(f"{args.output}_confusion.csv")
with open(f"{args.output}_metrics.txt", "w") as f:
    f.write(f"Adjusted Rand Index: {ari:.4f}\n")

print(f"✅ Comparison complete. ARI={ari:.4f}, outputs saved to {args.output}_*")
