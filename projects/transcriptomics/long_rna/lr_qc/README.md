# 🧾 Long RNA-seq Quality Control (QC) and Summary Assessment

## Overview
This module performs **post-analysis quality control** for long-read RNA-seq pipelines.  
It ensures that the **alignment**, **normalization**, and **differential expression** stages were performed correctly by checking sample quality, consistency, and reproducibility.

QC helps identify potential **outliers, batch effects, or low-quality samples** before final visualization and interpretation.  
It integrates statistical summaries, PCA-based inspection, and correlation heatmaps for robust assessment.

---

## 🔬 Pipeline Stage: Quality Control

**Purpose:**  
Validate the consistency and reliability of expression data after normalization or DE analysis.

**Key tasks:**
1. Evaluate expression distributions and gene detection rates.  
2. Visualize sample-level structure via PCA and clustering.  
3. Measure pairwise correlations between samples.  
4. Generate a unified QC report summarizing all visual outputs.

---

## ⚙️ Scripts and Utilities

| Script | Description | Input | Output | When to Use |
|--------|--------------|--------|---------|--------------|
| **summarize_expression_metrics.R** | Summarizes gene-level expression (counts per gene, detection rate). | Normalized counts matrix (TSV/CSV) | `gene_expression_summary.csv`, `gene_count_distribution.pdf` | Use right after normalization to check gene coverage and count distribution. |
| **run_pca_qc.R** | Performs PCA on normalized counts and plots sample clustering. | Normalized counts + metadata | `pca_samples.pdf` | Use to detect outliers or batch effects between biological replicates. |
| **run_correlation_heatmap.R** | Generates Pearson correlation heatmap between samples. | Normalized counts | `sample_correlation_heatmap.pdf` | Use to assess reproducibility across samples or replicates. |
| **generate_qc_report.R** | Combines all QC plots into a unified PDF report. | Directory containing QC PDFs | `qc_summary_report.pdf` | Use as the final summary step to compile QC visuals and statistics. |

---


## 📁 Outputs

| File | Description |
|------|--------------|
| gene_expression_summary.csv | Per-gene total counts and number of detected samples |
| gene_count_distribution.pdf | Histogram of total gene counts (log10 scale) |
| pca_samples.pdf | PCA plot showing sample clustering and condition grouping |
| sample_correlation_heatmap.pdf | Heatmap of pairwise sample correlations |
| qc_summary_report.pdf | Unified report compiling all QC plots |

---

## 🧠 Example Commands

Example 1: Summarize expression metrics  
Rscript summarize_expression_metrics.R --counts results/normalization/counts_tpm.csv --outdir results/qc

Example 2: PCA and sample clustering  
Rscript run_pca_qc.R --counts results/normalization/counts_tpm.csv --metadata data/sample_metadata.tsv --outdir results/qc

Example 3: Correlation heatmap  
Rscript run_correlation_heatmap.R --counts results/normalization/counts_tpm.csv --outdir results/qc

Example 4: Generate QC summary report  
Rscript generate_qc_report.R --indir results/qc --output results/qc_summary/qc_summary_report.pdf

---

