#!/usr/bin/env python3
"""
Script: run_preprocess_scanpy.py
Purpose: Preprocess scRNA-seq datasets with Scanpy (normalization, variable genes, scaling)
"""

import scanpy as sc
import argparse
import os

parser = argparse.ArgumentParser(description="Preprocess scRNA-seq with Scanpy")
parser.add_argument("-i", "--input", required=True, help="Comma-separated h5ad files")
parser.add_argument("-o", "--output", required=True, help="Output directory")
args = parser.parse_args()

inputs = args.input.split(",")
os.makedirs(args.output, exist_ok=True)

for file in inputs:
    print(f"Processing {file}")
    adata = sc.read_h5ad(file)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata)
    sc.pp.scale(adata)
    adata.write_h5ad(os.path.join(args.output, os.path.basename(file).replace(".h5ad","_processed.h5ad")))

print("Preprocessing finished for all samples.")
