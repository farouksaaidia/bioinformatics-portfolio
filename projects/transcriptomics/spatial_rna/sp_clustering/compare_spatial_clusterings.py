#!/usr/bin/env python3
"""
Compare spatial clustering results using Adjusted Rand Index and confusion matrices.
"""
import pandas as pd
import argparse
from sklearn.metrics import adjusted_rand_score, confusion_matrix

parser = argparse.ArgumentParser()
parser.add_argument("--a", required=True, help="CSV 1 with columns: spot_id,cluster")
parser.add_argument("--b", required=True, help="CSV 2 with columns: spot_id,cluster")
parser.add_argument("--out_prefix", required=True, help="Output prefix for comparison files")
args = parser.parse_args()

a = pd.read_csv(args.a)
b = pd.read_csv(args.b)
merged = pd.merge(a, b, on="spot_id", suffixes=("_a","_b"))

ari = adjusted_rand_score(merged.cluster_a, merged.cluster_b)
cm = pd.crosstab(merged.cluster_a, merged.cluster_b)

cm.to_csv(f"{args.out_prefix}_confusion.csv")
with open(f"{args.out_prefix}_metrics.txt","w") as f:
    f.write(f"Adjusted Rand Index (ARI): {ari:.4f}\n")

print(f"✅ Comparison done. ARI={ari:.4f}, results saved with prefix {args.out_prefix}")
