#!/usr/bin/env python3
"""
Map single-cell reference to spatial transcriptomics using Tangram.
"""
import scanpy as sc, tangram as tg
import argparse, os, sys

parser = argparse.ArgumentParser(description="Tangram-based spatial cell type mapping")
parser.add_argument("-r","--reference", required=True, help="Reference scRNA .h5ad")
parser.add_argument("-s","--spatial", required=True, help="Spatial .h5ad")
parser.add_argument("-o","--output", required=True, help="Output annotated spatial .h5ad")
args = parser.parse_args()

if not os.path.exists(args.reference) or not os.path.exists(args.spatial):
    sys.exit("❌ Missing input file(s)")

print("📥 Loading data...")
adata_sc = sc.read_h5ad(args.reference)
adata_sp = sc.read_h5ad(args.spatial)

print("🧠 Running Tangram mapping...")
tg.pp_adatas(adata_sc, adata_sp, genes=None)
ad_map = tg.map_cells_to_space(adata_sc, adata_sp, mode='clusters')
tg.project_cell_annotations(ad_map, adata_sp, annotation='cell_type')

adata_sp.obs['tangram_predicted'] = adata_sp.obs['cell_type_predicted']

print("💾 Saving annotated spatial object...")
adata_sp.write_h5ad(args.output)
print("✅ Tangram mapping complete ->", args.output)
