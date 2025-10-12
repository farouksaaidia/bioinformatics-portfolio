#!/usr/bin/env python3
import pandas as pd, argparse, os, sys
from itertools import combinations

parser = argparse.ArgumentParser(description="Compare DEG sets between clusters/conditions.")
parser.add_argument("-d","--deg_dir", required=True, help="Directory containing DEG CSV files")
parser.add_argument("-o","--output", required=True, help="Output directory")
args = parser.parse_args()
os.makedirs(args.output, exist_ok=True)

files = [f for f in os.listdir(args.deg_dir) if f.endswith(".csv")]
data = {os.path.splitext(f)[0]: pd.read_csv(os.path.join(args.deg_dir, f)) for f in files}

summary = []
for (n1, n2) in combinations(data.keys(), 2):
    g1 = set(data[n1]["gene"])
    g2 = set(data[n2]["gene"])
    overlap = len(g1 & g2)
    union = len(g1 | g2)
    jaccard = overlap / union if union > 0 else 0
    summary.append({"set1": n1, "set2": n2, "overlap": overlap, "jaccard": jaccard})
pd.DataFrame(summary).to_csv(os.path.join(args.output, "deg_overlap_summary.csv"), index=False)
print("✅ DEG set comparison complete.")
