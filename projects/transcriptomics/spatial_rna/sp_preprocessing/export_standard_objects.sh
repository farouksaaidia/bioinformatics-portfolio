#!/usr/bin/env bash
set -euo pipefail
# Exports canonical Seurat .rds to interop formats (.h5seurat -> .h5ad) for Python workflows.
# Usage: export_standard_objects.sh -r input_seurat.rds -o output_prefix

usage(){
  echo "Usage: $0 -r <input_seurat.rds> -o <output_prefix>"
  exit 1
}

INPUT_RDS=""
OUT_PREFIX=""

while getopts "r:o:" opt; do
  case $opt in
    r) INPUT_RDS=$OPTARG ;;
    o) OUT_PREFIX=$OPTARG ;;
    *) usage ;;
  esac
done

if [[ -z "$INPUT_RDS" || -z "$OUT_PREFIX" ]]; then
  usage
fi

if [[ ! -f "$INPUT_RDS" ]]; then
  echo "❌ Input file not found: $INPUT_RDS"
  exit 1
fi

R_SCRIPT=$(cat <<'RS'
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
})
args <- commandArgs(trailingOnly = TRUE)
input_rds <- args[1]
out_prefix <- args[2]
so <- readRDS(input_rds)
h5seurat <- paste0(out_prefix, '.h5seurat')
h5ad <- paste0(out_prefix, '.h5ad')
message('💾 Saving h5seurat -> ', h5seurat)
SaveH5Seurat(so, filename = h5seurat, overwrite = TRUE)
message('🔁 Converting to h5ad -> ', h5ad)
Convert(h5seurat, dest = "h5ad", overwrite = TRUE)
cat('✅ Export finished: ', h5seurat, ' and ', h5ad, '\n')
RS
)

Rscript -e "$R_SCRIPT" -- "$INPUT_RDS" "$OUT_PREFIX"
echo "✅ export_standard_objects.sh completed"
