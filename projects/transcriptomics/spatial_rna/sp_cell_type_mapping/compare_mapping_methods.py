#!/usr/bin/env python3
"""
Compare multiple spatial mapping methods (Seurat, Tangram, cell2location).
"""
import pandas as pd
import argparse
from sklearn.metrics import adjusted_rand_score

parser = argparse.ArgumentParser()
parser.add_argument("--a", required=True, help="CSV from method 1 (spot_id,cell_type)")
parser.add_argument("--b", required=True, help="CSV from method 2 (spot_id,cell_type)")
parser.add_argument("--out_prefix", required=True, help="Output file prefix")
args = parser.parse_args()

df1 = pd.read_csv(args.a)
df2 = pd.read_csv(args.b)
merged = pd.merge(df1, df2, on="spot_id", suffixes=("_a","_b"))

ari = adjusted_rand_score(merged.cell_type_a, merged.cell_type_b)
confmat = pd.crosstab(merged.cell_type_a, merged.cell_type_b)
confmat.to_csv(f"{args.out_prefix}_confusion.csv")

with open(f"{args.out_prefix}_metrics.txt", "w") as f:
    f.write(f"Adjusted Rand Index: {ari:.4f}\n")

print(f"✅ Comparison complete. ARI={ari:.4f}")
