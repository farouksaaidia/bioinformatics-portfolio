# sc_differential_expression

Differential expression and pathway analysis module for single-cell RNA-seq data.  
This stage follows **cell-type annotation** and identifies biologically relevant genes and pathways distinguishing clusters, conditions, or states.

---

## 📂 Folder Purpose

`projects/transcriptomics/sc_rna/sc_differential_expression/`  
Contains all scripts related to **differential gene expression**, **functional enrichment**, **comparisons**, and **reporting**.

---

## 🧩 Scripts Overview

| Script | Language | Input | Output | Description |
|--------|-----------|--------|---------|-------------|
| **find_differential_genes.R** | R | annotated Seurat `.rds` | per-cluster DEG CSVs | Finds DEGs per cluster or condition using Seurat `FindMarkers()` |
| **visualize_deg_results.R** | R | Seurat `.rds`, DEG CSVs | volcano, heatmap images | Generates volcano and heatmap plots of significant genes |
| **pathway_enrichment_analysis.R** | R | DEG CSVs | GO/KEGG enrichment CSVs | Performs GO and KEGG enrichment via `clusterProfiler` |
| **compare_deg_sets.py** | Python | DEG CSVs | overlap table | Quantifies overlap (Jaccard index) between DEG sets |
| **generate_deg_report.R** | R | DEG & enrichment results | summary PDF | Produces combined visual PDF report of differential analysis |
| **.gitkeep** | text | — | — | Keeps folder tracked when empty |

---

## ⚙️ Usage Examples

### 1️⃣ Differential Gene Detection
```bash
projects/transcriptomics/sc_rna/sc_differential_expression/find_differential_genes.R \
  -i results/annotated_seurat.rds \
  -g cell_type \
  -o results/differential_expression/DEG_tables
```

### 2️⃣ Visualize DEGs
```bash
projects/transcriptomics/sc_rna/sc_differential_expression/visualize_deg_results.R \
  -i results/annotated_seurat.rds \
  -d results/differential_expression/DEG_tables \
  -o results/differential_expression/plots
```

### 3️⃣ Functional Enrichment
```bash
projects/transcriptomics/sc_rna/sc_differential_expression/pathway_enrichment_analysis.R \
  -d results/differential_expression/DEG_tables \
  -o results/differential_expression/enrichment
```

### 4️⃣ Compare DEG Sets
```bash
projects/transcriptomics/sc_rna/sc_differential_expression/compare_deg_sets.py \
  -d results/differential_expression/DEG_tables \
  -o results/differential_expression/comparisons
```

### 5️⃣ Generate Summary Report
```bash
projects/transcriptomics/sc_rna/sc_differential_expression/generate_deg_report.R \
  -d results/differential_expression/DEG_tables \
  -e results/differential_expression/enrichment \
  -o reports/differential_expression_summary.pdf
```

---

## 🧠 Best Practices

- Always use **normalized and scaled** data for DEG analysis.  
- Ensure grouping metadata (e.g., `cell_type`, `condition`) is accurate and non-empty.  
- For multi-condition experiments, consider including batch-corrected data or using models (e.g., MAST, DESeq2) for higher precision.  
- Apply multiple testing correction (`p_val_adj`) and consider only DEGs with `|log2FC| > 0.25` and `FDR < 0.05`.  
- Use biologically relevant **backgrounds** when running enrichment (e.g., expressed genes).  
- Validate enriched pathways using multiple ontology sources (GO, KEGG, Reactome).  
- Record versions of all R packages (`sessionInfo()`) for reproducibility.  
- Prefer interactive exploration of DEGs in later visualization notebooks.

---

## 🧾 Attribution
Created and maintained by **Farouk Saaidia (2025)**.  
For research use, cite as:

> Saaidia F. (2025). *Differential Expression & Pathway Analysis Module for Single-Cell RNA-Seq.*

