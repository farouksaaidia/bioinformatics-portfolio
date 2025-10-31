#!/usr/bin/env python3
import pandas as pd
import argparse
import sys

parser = argparse.ArgumentParser(description="Compare DE gene lists across DESeq2, edgeR, and IsoformSwitchAnalyzeR results.")
parser.add_argument("-d", "--deseq2", required=True, help="DESeq2 results CSV")
parser.add_argument("-e", "--edger", required=True, help="edgeR results CSV")
parser.add_argument("-i", "--isoform", required=False, help="IsoformSwitchAnalyzeR results CSV")
parser.add_argument("-o", "--output", required=True, help="Output summary TSV")
args = parser.parse_args()

try:
    deseq = pd.read_csv(args.deseq2, index_col=0)
    edger = pd.read_csv(args.edger, index_col=0)
    isoform = pd.read_csv(args.isoform, index_col=0) if args.isoform else None
except Exception as e:
    sys.exit(f"❌ Error loading input files: {e}")

summary = pd.DataFrame({
    "DESeq2_sig": (deseq["padj"] < 0.05).sum(),
    "edgeR_sig": (edger["FDR"] < 0.05).sum(),
    "Overlap": len(set(deseq.index[deseq["padj"] < 0.05]) & set(edger.index[edger["FDR"] < 0.05]))
}, index=[0])

if isoform is not None:
    summary["IsoformSwitch_sig"] = (isoform["isoform_switch_q_value"] < 0.05).sum()

summary.to_csv(args.output, sep="\t", index=False)
print("✅ Comparison complete. Results saved to", args.output)
