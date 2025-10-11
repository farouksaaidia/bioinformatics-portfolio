#!/usr/bin/env python3
"""
Ensemble annotation:
- Inputs: scanpy .h5ad annotated by CellTypist (obs key: predicted_labels), a CSV marker-based annotation (cell_id, cell_type),
  and optionally a SingleR-style CSV (cell_id, cell_type_singleR).
- Outputs: .h5ad with ensemble label in obs['ensemble_label'] and a CSV with per-cell votes/confidence.
Consensus rules:
- If >=2 methods agree -> that label
- If all different -> label 'ambiguous'
- Produces confidence score = fraction of methods agreeing
"""
import argparse
import os
import sys
import pandas as pd
import scanpy as sc

parser = argparse.ArgumentParser(description="Ensemble annotation aggregator")
parser.add_argument("-i","--input_h5ad", required=True, help="Input .h5ad with at least CellTypist obs key or predictions")
parser.add_argument("--celltypist_key", default="predicted_labels", help="obs key for CellTypist predictions (default: predicted_labels)")
parser.add_argument("--marker_csv", default=None, help="CSV from marker-based method (columns: cell_id, cell_type)")
parser.add_argument("--singler_csv", default=None, help="CSV from SingleR method (columns: cell_id, cell_type_singleR)")
parser.add_argument("-o","--output_h5ad", required=True, help="Output annotated .h5ad")
parser.add_argument("--out_votes", default=None, help="Optional CSV to write per-cell votes")
args = parser.parse_args()

if not os.path.exists(args.input_h5ad):
    sys.exit(f"❌ {args.input_h5ad} not found")

adata = sc.read_h5ad(args.input_h5ad)
cells = adata.obs_names.to_list()
votes = pd.DataFrame(index=cells)

methods = []
# CellTypist
if args.celltypist_key in adata.obs.columns:
    votes['celltypist'] = adata.obs[args.celltypist_key].astype(str).values
    methods.append('celltypist')
else:
    votes['celltypist'] = None

# Marker-based
if args.marker_csv:
    if not os.path.exists(args.marker_csv):
        sys.exit(f"❌ marker_csv {args.marker_csv} not found")
    mdf = pd.read_csv(args.marker_csv)
    if not {'cell_id','cell_type'}.issubset(mdf.columns):
        sys.exit("❌ marker_csv must have columns: cell_id, cell_type")
    mdf = mdf.set_index('cell_id').reindex(cells)
    votes['marker'] = mdf['cell_type'].astype(str).values
    methods.append('marker')

# SingleR
if args.singler_csv:
    if not os.path.exists(args.singler_csv):
        sys.exit(f"❌ singler_csv {args.singler_csv} not found")
    sdf = pd.read_csv(args.singler_csv)
    if not {'cell_id','cell_type_singleR'}.issubset(sdf.columns):
        sys.exit("❌ singler_csv must have columns: cell_id, cell_type_singleR")
    sdf = sdf.set_index('cell_id').reindex(cells)
    votes['singler'] = sdf['cell_type_singleR'].astype(str).values
    methods.append('singler')

if len(methods) == 0:
    sys.exit("❌ No annotation methods provided. Provide at least celltypist results inside the .h5ad or marker/singler CSVs")

# compute consensus
def consensus_label(row):
    vals = [v for v in row.values if (v is not None and v != 'nan' and v != 'None')]
    if len(vals)==0:
        return pd.Series({"label":"unknown","confidence":0.0})
    counts = pd.Series(vals).value_counts()
    top = counts.index[0]
    conf = counts.iloc[0] / len(vals)
    if conf >= 0.5 and counts.iloc[0] >= 2:
        return pd.Series({"label":top,"confidence":conf})
    elif conf == 1.0:
        return pd.Series({"label":top,"confidence":conf})
    else:
        return pd.Series({"label":"ambiguous","confidence":conf})

cons = votes.apply(consensus_label, axis=1)
votes['ensemble_label'] = cons['label']
votes['ensemble_confidence'] = cons['confidence']

# write back to adata
adata.obs['ensemble_label'] = votes['ensemble_label'].values
adata.obs['ensemble_confidence'] = votes['ensemble_confidence'].values

adata.write_h5ad(args.output_h5ad)
if args.out_votes:
    votes.reset_index().rename(columns={'index':'cell_id'}).to_csv(args.out_votes, index=False)

print(f"✅ Ensemble annotation saved to {args.output_h5ad}")
if args.out_votes:
    print(f"✅ Votes saved to {args.out_votes}")
