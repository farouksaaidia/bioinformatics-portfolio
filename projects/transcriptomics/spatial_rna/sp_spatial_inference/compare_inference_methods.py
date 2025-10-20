#!/usr/bin/env python3
"""
Compare ligand–receptor inference results from multiple methods.
"""
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--nichenet", required=False, help="NicheNet output CSV")
parser.add_argument("--liana", required=False, help="LIANA output CSV")
parser.add_argument("--giotto", required=False, help="Giotto neighborhood output CSV")
parser.add_argument("--out_prefix", required=True, help="Output file prefix")
args = parser.parse_args()

dfs = []
if args.nichenet: dfs.append(pd.read_csv(args.nichenet).assign(method="NicheNet"))
if args.liana: dfs.append(pd.read_csv(args.liana).assign(method="LIANA"))
if args.giotto: dfs.append(pd.read_csv(args.giotto).assign(method="Giotto"))

if not dfs:
    raise SystemExit("❌ No input files provided")

merged = pd.concat(dfs, ignore_index=True)
common_pairs = merged.groupby(["ligand","receptor"]).method.nunique().reset_index()
common_pairs = common_pairs[common_pairs.method > 1]

merged.to_csv(f"{args.out_prefix}_merged.csv", index=False)
common_pairs.to_csv(f"{args.out_prefix}_common.csv", index=False)
print("✅ Comparison complete. Results saved with prefix", args.out_prefix)
