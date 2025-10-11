#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 -r <input_rds> -h <input_h5ad> -o <output_dir>"
  echo "Either -r or -h must be provided (or both)."
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

if [[ -z "${RDS}" && -z "${H5AD}" ]]; then
  echo "❌ Provide at least -r or -h"
  usage
fi

if [[ -z "${OUTDIR}" ]]; then
  echo "❌ Provide -o <output_dir>"
  usage
fi

mkdir -p "${OUTDIR}"

manifest="${OUTDIR}/export_manifest.tsv"
echo -e "source\tfilename\tsize_bytes" > "${manifest}"

if [[ -n "${RDS}" && -f "${RDS}" ]]; then
  cp "${RDS}" "${OUTDIR}/annotated_seurat.rds"
  sz=$(stat -c%s "${OUTDIR}/annotated_seurat.rds")
  echo -e "${RDS}\tannotated_seurat.rds\t${sz}" >> "${manifest}"
  echo "✅ Exported annotated Seurat object to ${OUTDIR}/annotated_seurat.rds"
elif [[ -n "${RDS}" ]]; then
  echo "⚠️ RDS file not found: ${RDS}"
fi

if [[ -n "${H5AD}" && -f "${H5AD}" ]]; then
  cp "${H5AD}" "${OUTDIR}/annotated_scanpy.h5ad"
  sz=$(stat -c%s "${OUTDIR}/annotated_scanpy.h5ad")
  echo -e "${H5AD}\tannotated_scanpy.h5ad\t${sz}" >> "${manifest}"
  echo "✅ Exported annotated Scanpy object to ${OUTDIR}/annotated_scanpy.h5ad"
elif [[ -n "${H5AD}" ]]; then
  echo "⚠️ H5AD file not found: ${H5AD}"
fi

echo "📄 Manifest saved to ${manifest}"
echo "✅ Export complete."
