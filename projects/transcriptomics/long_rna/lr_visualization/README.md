# 📊 Long RNA-seq Visualization and Reporting

## Overview
This module generates **publication-ready visualizations** and **interactive dashboards** summarizing the outcomes of normalization, differential expression, and enrichment analyses.  
It provides a cohesive visual interpretation of the data—highlighting significantly altered genes, pathways, and biological mechanisms revealed in earlier stages.

The visualization layer ensures that key findings are communicated clearly and that QC, DE, and enrichment results are intuitively interpretable.

---

## 🔬 Pipeline Stage: Visualization and Reporting

**Purpose:**  
Transform statistical results into high-quality visual representations suitable for reports, presentations, and publications.

**Key tasks:**
1. Visualize differential expression through volcano and heatmap plots.  
2. Summarize enriched biological terms and pathways in bubble plots.  
3. Provide an interactive Shiny dashboard for real-time exploration.  
4. Compile a unified report integrating all generated plots.

---

## ⚙️ Scripts and Utilities

| Script | Description | Input | Output | When to Use |
|--------|--------------|--------|---------|--------------|
| **plot_volcano_de.R** | Generates a volcano plot highlighting up/down-regulated genes. | DE results (`gene_id`, `log2FoldChange`, `padj`) | `volcano_plot.pdf` | Use to visualize global differential expression and identify most significant genes. |
| **plot_heatmap_topgenes.R** | Draws a heatmap of top DE genes based on adjusted p-value. | Normalized counts, DE results | `heatmap_top_genes.pdf` | Use to visualize clustering patterns of the top N differentially expressed genes. |
| **plot_enrichment_bubble.R** | Plots top enriched terms or pathways using bubble visualization. | Enrichment results (GO/KEGG/Reactome) | `enrichment_bubble_plot.pdf` | Use to visually summarize functional enrichment outcomes. |
| **interactive_dashboard.R** | Launches a Shiny app for interactive exploration of DE results. | DE result CSV | Local interactive dashboard | Use to explore fold changes and significance interactively. |
| **generate_visual_report.R** | Combines all static plots (PDFs) into a single visual report. | Directory with visualization outputs | `visualization_summary_report.pdf` | Use as the final reporting step to compile all plots into one file. |

---

## 🧠 Example Commands

Example 1: Volcano plot  
Rscript plot_volcano_de.R --input results/deseq2/deseq2_results.csv --outdir results/visualization

Example 2: Heatmap of top DE genes  
Rscript plot_heatmap_topgenes.R --counts results/normalization/counts_tpm.csv --de results/deseq2/deseq2_results.csv --outdir results/visualization --topn 50

Example 3: Enrichment bubble plot  
Rscript plot_enrichment_bubble.R --input results/enrichment/go_enrichment_results.csv --outdir results/visualization

Example 4: Interactive dashboard  
Rscript interactive_dashboard.R  
(then upload DE results CSV via browser interface)

Example 5: Generate visual summary report  
Rscript generate_visual_report.R --indir results/visualization --output results/final_reports/visualization_summary.pdf

---

## 📁 Outputs

| File | Description |
|------|--------------|
| volcano_plot.pdf | Volcano plot of differential expression |
| heatmap_top_genes.pdf | Heatmap of top DE genes |
| enrichment_bubble_plot.pdf | Bubble chart of top enriched terms |
| visualization_summary_report.pdf | Combined visual report of all generated plots |
| Interactive Dashboard | Shiny app interface for dynamic exploration of DE data |

---

## 🧩 Dependencies
- **R packages:** optparse, ggplot2, pheatmap, shiny, readr, dplyr, gridExtra, ggrepel, png  
- **Runtime:** R ≥ 4.0.0  
- **Input requirement:** Differential expression and enrichment result files.

---


This stage consolidates the entire analysis into an accessible, visual format suitable for interpretation, collaboration, and publication.

