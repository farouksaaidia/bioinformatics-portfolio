#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<USAGE
Usage: $0 -r <input_rds> -h <input_h5ad> -o <output_dir>

Copies annotated Seurat (.rds) and Scanpy (.h5ad) objects into an output directory.
Examples:
  $0 -r results/annotated_seurat.rds -h results/annotated_scanpy.h5ad -o exports/
USAGE
  exit 1
}

RDS=""
H5AD=""
OUTDIR=""

while getopts "r:h:o:" opt; do
  case $opt in
    r) RDS=$OPTARG ;;
    h) H5AD=$OPTARG ;;
    o) OUTDIR=$OPTARG ;;
    *) usage ;;
  esac
done

if [[ -z "${OUTDIR}" ]]; then
  echo "❌ Missing -o <output_dir>"
  usage
fi

mkdir -p "$OUTDIR"

if [[ -n "${RDS}" && -f "${RDS}" ]]; then
  cp "$RDS" "$OUTDIR/annotated_seurat.rds"
  echo "✅ Exported annotated Seurat object to $OUTDIR/annotated_seurat.rds"
else
  echo "⚠️ Seurat RDS not provided or not found: ${RDS}"
fi

if [[ -n "${H5AD}" && -f "${H5AD}" ]]; then
  cp "$H5AD" "$OUTDIR/annotated_scanpy.h5ad"
  echo "✅ Exported annotated Scanpy object to $OUTDIR/annotated_scanpy.h5ad"
else
  echo "⚠️ Scanpy H5AD not provided or not found: ${H5AD}"
fi

echo "✅ Export step finished."
