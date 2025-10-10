#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 -r <input_rds> -h <input_h5ad> -o <output_dir>"
  exit 1
}

while getopts "r:h:o:" opt; do
  case $opt in
    r) RDS=$OPTARG ;;
    h) H5AD=$OPTARG ;;
    o) OUTDIR=$OPTARG ;;
    *) usage ;;
  esac
done

mkdir -p "$OUTDIR"

if [[ -f "$RDS" ]]; then
  cp "$RDS" "$OUTDIR/annotated_seurat.rds"
  echo "✅ Exported annotated Seurat object to $OUTDIR/annotated_seurat.rds"
fi

if [[ -f "$H5AD" ]]; then
  cp "$H5AD" "$OUTDIR/annotated_scanpy.h5ad"
  echo "✅ Exported annotated Scanpy object to $OUTDIR/annotated_scanpy.h5ad"
fi
