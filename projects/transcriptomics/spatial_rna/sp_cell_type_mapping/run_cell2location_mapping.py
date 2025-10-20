#!/usr/bin/env python3
"""
Probabilistic cell type mapping using cell2location.
"""
import scvi, cell2location, scanpy as sc
import argparse, os, sys

parser = argparse.ArgumentParser()
parser.add_argument("-r","--reference", required=True, help="Reference .h5ad with cell types")
parser.add_argument("-s","--spatial", required=True, help="Spatial .h5ad")
parser.add_argument("-o","--output", required=True, help="Output annotated .h5ad")
args = parser.parse_args()

print("📥 Loading datasets...")
adata_sc = sc.read_h5ad(args.reference)
adata_sp = sc.read_h5ad(args.spatial)

print("⚙️ Training reference model...")
scvi.model.SCVI.setup_anndata(adata_sc, labels_key="cell_type")
model = scvi.model.SCVI(adata_sc)
model.train()

print("🔗 Running cell2location mapping...")
cell2location.models.Cell2location.setup_anndata(adata_sp)
c2l_model = cell2location.models.Cell2location.from_scvi_model(model, adata_sp)
c2l_model.train(max_epochs=200)
adata_sp = c2l_model.export_posterior(adata_sp)

adata_sp.write_h5ad(args.output)
print("✅ Cell2location mapping complete ->", args.output)
