# 🎯 Long RNA-seq Differential Expression Analysis

## Overview
This module identifies significantly differentially expressed **genes** and **isoforms** across experimental conditions using multiple statistical frameworks.  
It supports **count-based** (DESeq2, edgeR) and **isoform-level** (IsoformSwitchAnalyzeR) analyses, ensuring robust inference for both gene-level and transcript-level dynamics.

Differential expression (DE) analysis compares expression values between groups (e.g., treatment vs. control) to detect transcripts or genes that change significantly under different conditions.  
The results can be used for pathway enrichment, biomarker discovery, and biological interpretation of transcriptional changes.

---

## 🔬 Pipeline Stage: Differential Expression

**Purpose:**  
Identify transcripts or genes with statistically significant expression differences between sample groups.

**Key tasks:**
1. Statistical modeling and normalization-aware DE testing (DESeq2 / edgeR).  
2. Isoform-level DE testing and switch detection (IsoformSwitchAnalyzeR).  
3. Comparison of DE results between methods for validation.  
4. Generation of ranked DE tables for downstream functional enrichment.

---

## ⚙️ Scripts and Utilities

| Script | Description | Input | Output | When to Use |
|--------|--------------|--------|---------|--------------|
| **run_deseq2.R** | Performs DE analysis using DESeq2 (negative binomial model). | Normalized count matrix (CSV/TSV), metadata (TSV with sample-condition) | `deseq2_results.csv`, `normalized_counts.csv` | Use when you want robust DE at the **gene level** with moderate sample sizes (≥3 per group). |
| **run_edger.R** | Runs DE analysis using edgeR with quasi-likelihood tests. | Raw counts (TSV), metadata (TSV: sample, condition) | `edger_results.csv` | Use for large cohorts or complex designs; performs well with low-count genes. |
| **run_isoform_switch_analysis.R** | Detects isoform-level DE and isoform switching events using IsoformSwitchAnalyzeR. | Isoform abundances (TPM/Counts), metadata, annotation GTF | `isoform_switch_summary.csv` | Use for long-read datasets or when isoform structure and alternative splicing are biologically important. |
| **compare_de_results.py** | Compares DESeq2, edgeR, and IsoformSwitch results; outputs overlap metrics. | DE result tables (CSV/TSV) | `de_comparison_summary.tsv` | Use to validate consistency across methods and find commonly regulated genes. |

---

## 🧠 Example Commands

Example 1: DESeq2 analysis
Rscript run_deseq2.R --counts data/normalized_counts.csv --metadata data/sample_metadata.tsv --outdir results/deseq2

Example 2: edgeR analysis
Rscript run_edger.R --counts data/gene_counts.tsv --metadata data/sample_metadata.tsv --outdir results/edger

Example 3: Isoform-level DE
Rscript run_isoform_switch_analysis.R --abundance data/isoform_abundance.tsv --metadata data/sample_metadata.tsv --annotation data/annotation.gtf --outdir results/isoform_switch

Example 4: Compare DE results across methods
python3 compare_de_results.py --deseq2 results/deseq2/deseq2_results.csv --edger results/edger/edger_results.csv --isoform results/isoform_switch/isoform_switch_summary.csv --output results/de_comparison_summary.tsv

---

## 📁 Outputs

| File | Description |
|------|--------------|
| deseq2_results.csv | Differential expression results from DESeq2 |
| edger_results.csv | Differential expression results from edgeR |
| isoform_switch_summary.csv | Isoform-level switch and DE summary |
| de_comparison_summary.tsv | Overlap and metrics between DE methods |
| normalized_counts.csv | Normalized count matrix from DESeq2 |

---

## 🧩 Dependencies
- **R packages:** DESeq2, edgeR, IsoformSwitchAnalyzeR, optparse, readr  
- **Python packages:** pandas  
- **Runtime:** R ≥ 4.0.0, Python ≥ 3.8

---

## 🧭 Next Step
Proceed to [`lr_enrichment`](../lr_enrichment/), where DE gene or isoform lists are used for functional enrichment analysis (GO, KEGG, Reactome).

