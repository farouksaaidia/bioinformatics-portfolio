#!/usr/bin/env python3
import argparse, os, sys, pandas as pd, glob

parser = argparse.ArgumentParser(description="Merge Salmon quant.sf files into count and TPM matrices")
parser.add_argument("-i","--indir", required=True, help="Directory containing sample quant folders")
parser.add_argument("-o","--outdir", required=True, help="Output directory")
args = parser.parse_args()

if not os.path.isdir(args.indir):
    sys.exit("Input directory not found")

samples = sorted([d for d in os.listdir(args.indir) if os.path.isdir(os.path.join(args.indir, d))])
if len(samples) == 0:
    sys.exit("No Salmon quant directories found")

os.makedirs(args.outdir, exist_ok=True)

count_dfs = []
tpm_dfs   = []

for s in samples:
    quant = os.path.join(args.indir, s, "quant.sf")
    if not os.path.exists(quant):
        sys.exit(f"quant.sf missing for {s}")
    df = pd.read_csv(quant, sep="\t")
    count_dfs.append(df[["Name","NumReads"]].rename(columns={"NumReads":s}))
    tpm_dfs.append(df[["Name","TPM"]].rename(columns={"TPM":s}))

counts = count_dfs[0]
for df in count_dfs[1:]:
    counts = counts.merge(df, on="Name")

tpms = tpm_dfs[0]
for df in tpm_dfs[1:]:
    tpms = tpms.merge(df, on="Name")

counts.to_csv(os.path.join(args.outdir, "salmon_counts_matrix.csv"), index=False)
tpms.to_csv(os.path.join(args.outdir, "salmon_tpm_matrix.csv"), index=False)

print("Merged Salmon matrices saved:")
print(" - salmon_counts_matrix.csv")
print(" - salmon_tpm_matrix.csv")
