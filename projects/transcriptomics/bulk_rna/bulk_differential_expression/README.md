# Bulk RNA-seq Visualization

This folder contains scripts to generate different visualizations for bulk RNA-seq results.  
These visualizations help interpret sample quality, variability, and biological significance of differential expression.

---

## Visualization Scripts

| Script | Input | Output | When to Use | Particularity |
|--------|-------|--------|-------------|---------------|
| `run_pca_plot.R` | Normalized count matrix, sample metadata | `pca_plot.png` | To check global sample clustering and batch effects | Highlights variation across conditions |
| `run_heatmap.R` | Normalized count matrix | `heatmap_topgenes.png` | To visualize top variable genes across samples | Shows gene expression similarity/differences |
| `run_ma_plot.R` | DE results (`baseMean`, `log2FC`, `padj`) | `ma_plot.png` | To assess systematic expression differences | Visualizes fold change vs expression strength |
| `run_volcano_plot.R` | DE results | `volcano_plot.png` | To quickly identify significant up/down-regulated genes | Combines fold change + statistical significance |
| `run_sample_distance_heatmap.R` | Normalized count matrix, metadata | `sample_distance_heatmap.png` | To evaluate similarity between samples | Clustering based on expression distance |
| `run_barplot_topgenes.R` | DE results | `barplot_topgenes.png` | To show most significant DEGs | Focus on top candidates |
| `run_pathway_enrichment_plot.R` | Pathway/GO enrichment results | `pathway_enrichment_plot.png` | To interpret biological functions of DEGs | Aggregates DEGs into pathways |
| `run_gsea_plot.R` | Ranked gene list, pathway GMT file | `gsea_plot.png` | To explore enriched gene sets without thresholds | Detects subtle but coordinated changes |
| `run_boxplot_genes.R` | Expression matrix, metadata, gene of interest | `boxplot_<gene>.png` | To validate expression of specific candidate genes | Gene-centric visualization |
| `run_upset_plot.R` | Binary DEG membership matrix | `upset_plot.png` | To compare overlap of DEGs across methods/conditions | Better than Venn for many comparisons |

---

## Notes
- Use **PCA, heatmaps, and distance plots** early to assess QC and batch effects.  
- Use **MA, volcano, and barplots** after DE analysis to highlight candidate genes.  
- Use **Pathway and GSEA plots** for biological interpretation.  
- **UpSet plots** are useful if multiple DEG methods or comparisons are run.  

