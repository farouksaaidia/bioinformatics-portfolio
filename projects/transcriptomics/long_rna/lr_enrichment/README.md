# 🧬 Long RNA-seq Functional Enrichment Analysis

## Overview
This module performs **functional enrichment analysis** on differentially expressed (DE) genes or isoforms.  
The goal is to interpret DE results biologically by identifying overrepresented **Gene Ontology (GO) terms**, **KEGG pathways**, and **Reactome pathways**, and to perform **GSEA (Gene Set Enrichment Analysis)** for ranked expression profiles.

Functional enrichment highlights biological processes, molecular functions, and pathways that are statistically enriched among DE genes, providing insights into the biological mechanisms underlying observed expression changes.

---

## 🔬 Pipeline Stage: Functional Enrichment

**Purpose:**  
Translate DE results into biological meaning by mapping significant genes to pathways and ontologies.

**Key tasks:**
1. Identify enriched GO terms (biological processes, molecular functions).  
2. Detect KEGG and Reactome pathway enrichment.  
3. Perform GSEA using ranked log2 fold changes.  
4. Generate visual summaries of top enriched terms and pathways.

---

## ⚙️ Scripts and Utilities

| Script | Description | Input | Output | When to Use |
|--------|--------------|--------|---------|--------------|
| **run_go_enrichment.R** | Performs Gene Ontology enrichment (BP, MF, CC) using clusterProfiler. | DE results (`gene_id`, `padj`) | `go_enrichment_results.csv` | Use to identify biological processes affected by differential expression. |
| **run_kegg_enrichment.R** | Runs KEGG pathway enrichment using Entrez/KEGG database. | DE results (`gene_id`, `padj`) | `kegg_enrichment_results.csv` | Use to explore pathway-level functional changes (metabolic, signaling). |
| **run_reactome_enrichment.R** | Performs Reactome pathway enrichment analysis. | DE results (`gene_id`, `padj`) | `reactome_enrichment_results.csv` | Use when Reactome-level annotations are preferred (human-focused, curated). |
| **run_gsea_analysis.R** | Executes ranked GSEA based on log2 fold changes using GO terms. | DE results (`gene_id`, `log2FoldChange`, `padj`) | `gsea_results.csv` | Use when you want to capture graded expression trends, not just thresholds. |
| **generate_enrichment_report.R** | Combines GO/KEGG/Reactome/GSEA results and visualizes top terms. | Enrichment result files | PDF plots of top terms/pathways | Use to summarize all enrichment outcomes into visual reports. |

---

## 🧠 Example Commands

Example 1: GO enrichment  
Rscript run_go_enrichment.R --input results/deseq2/deseq2_results.csv --outdir results/enrichment/go

Example 2: KEGG enrichment  
Rscript run_kegg_enrichment.R --input results/deseq2/deseq2_results.csv --outdir results/enrichment/kegg

Example 3: Reactome enrichment  
Rscript run_reactome_enrichment.R --input results/deseq2/deseq2_results.csv --outdir results/enrichment/reactome

Example 4: GSEA analysis  
Rscript run_gsea_analysis.R --input results/deseq2/deseq2_results.csv --outdir results/enrichment/gsea

Example 5: Generate global enrichment report  
Rscript generate_enrichment_report.R --indir results/enrichment/ --outdir results/enrichment_summary

---

## 📁 Outputs

| File | Description |
|------|--------------|
| go_enrichment_results.csv | GO Biological Process/Molecular Function enrichment results |
| kegg_enrichment_results.csv | KEGG pathway enrichment results |
| reactome_enrichment_results.csv | Reactome pathway enrichment results |
| gsea_results.csv | Gene Set Enrichment Analysis results |
| *_top10.pdf | Barplots of top enriched pathways/terms |
| enrichment_summary/ | Folder with combined enrichment visualizations |

---

## 🧩 Dependencies
- **R packages:** clusterProfiler, ReactomePA, org.Hs.eg.db, optparse, ggplot2, dplyr, readr  
- **Runtime:** R ≥ 4.0.0  
- **Input requirement:** Differential expression results containing at least `gene_id`, `log2FoldChange`, and `padj` columns.

---

## 🧭 Next Step
Proceed to [`lr_qc`](../lr_qc/), where multi-level quality control of normalized and enriched datasets is performed to ensure data integrity before final visualization and reporting.

