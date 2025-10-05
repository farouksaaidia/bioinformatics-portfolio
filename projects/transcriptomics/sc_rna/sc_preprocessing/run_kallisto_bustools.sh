#!/bin/bash
# =========================================
# run_kallisto_bustools.sh
# Preprocessing scRNA-seq using Kallisto + Bustools
# =========================================

set -euo pipefail

usage() {
    echo "Usage: $0 -i <sample_id> -1 <R1.fastq.gz> -2 <R2.fastq.gz> -i <index_file> -o <output_dir>"
    exit 1
}

while getopts ":s:1:2:i:o:h" opt; do
    case $opt in
        s) SAMPLE=$OPTARG ;;
        1) R1=$OPTARG ;;
        2) R2=$OPTARG ;;
        i) INDEX=$OPTARG ;;
        o) OUTDIR=$OPTARG ;;
        h) usage ;;
    esac
done

if [[ -z "${SAMPLE:-}" || -z "${R1:-}" || -z "${R2:-}" || -z "${INDEX:-}" || -z "${OUTDIR:-}" ]]; then
    echo "[ERROR] Missing arguments"; usage
fi

echo "[$(date)] Running Kallisto|Bustools for $SAMPLE"
kallisto bus -i "$INDEX" -x 10xv2 -o "$OUTDIR/$SAMPLE" "$R1" "$R2"
cd "$OUTDIR/$SAMPLE"
bustools correct -w whitelist.txt -o output.corrected.bus output.bus
bustools sort -o output.sorted.bus output.corrected.bus
bustools count -o counts -g transcripts_to_genes.txt -e matrix.ec -t transcripts.txt output.sorted.bus
echo "[$(date)] Finished Kallisto|Bustools for $SAMPLE"
