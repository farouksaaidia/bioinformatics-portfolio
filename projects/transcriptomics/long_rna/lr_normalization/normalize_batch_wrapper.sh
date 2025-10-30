#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 -c <counts.tsv> -l <lengths.tsv or empty> -M <manifest_or_single> -o <outdir> -m <methods>"
  echo "Manifest format (optional): a TSV with header sample[TAB]path_to_counts_column (if multiple sample-level files are used)."
  exit 1
}

COUNTS=""
LENGTHS=""
MANIFEST=""
OUTDIR=""
METHODS="tpm,cpm,tmm"

while getopts "c:l:M:o:m:" opt; do
  case $opt in
    c) COUNTS=$OPTARG ;;
    l) LENGTHS=$OPTARG ;;
    M) MANIFEST=$OPTARG ;;
    o) OUTDIR=$OPTARG ;;
    m) METHODS=$OPTARG ;;
    *) usage ;;
  esac
done

[[ -z "$COUNTS" || -z "$OUTDIR" ]] && usage
mkdir -p "$OUTDIR"

# If manifest provided, iterate (manifest expected to list counts files per sample or a single counts matrix)
if [[ -n "${MANIFEST}" && -f "${MANIFEST}" ]]; then
  echo "🔁 Running normalization per manifest entry (manifest: $MANIFEST)"
  while IFS=$'\t' read -r sample path; do
    echo "→ Processing sample: $sample (counts: $path)"
    sample_out="$OUTDIR/$sample"
    mkdir -p "$sample_out"
    Rscript projects/transcriptomics/long_rna/lr_normalization/compute_normalizations.R -i "$path" $( [[ -n "$LENGTHS" ]] && echo "-l $LENGTHS" ) -m "$METHODS" -o "$sample_out"
  done < <(tail -n +2 "$MANIFEST")
else
  echo "🔁 Running normalization on provided counts matrix"
  Rscript projects/transcriptomics/long_rna/lr_normalization/compute_normalizations.R -i "$COUNTS" $( [[ -n "$LENGTHS" ]] && echo "-l $LENGTHS" ) -m "$METHODS" -o "$OUTDIR"
fi

echo "✅ Normalization wrapper finished. Outputs in $OUTDIR"
