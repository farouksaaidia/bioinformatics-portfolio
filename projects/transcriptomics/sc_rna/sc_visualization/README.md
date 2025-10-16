# sc_visualization

Visualization and exploratory analysis module for single-cell RNA-seq data.  
This module provides both **publication-quality static plots** and **interactive dashboards** to explore annotated single-cell datasets across cell types, conditions, and genes.

---

## 📂 Folder Purpose

`projects/transcriptomics/sc_rna/sc_visualization/`  
Contains scripts for generating embeddings, feature plots, composition analyses, interactive dashboards, and summary reports.

---

## 🧩 Scripts Overview

| Script | Language | Input | Output | Description |
|--------|-----------|--------|---------|-------------|
| **plot_embeddings.R** | R | Seurat `.rds` | UMAP/tSNE plots (PNG/PDF) | Generates low-dimensional embeddings colored by metadata |
| **plot_features.R** | R | Seurat `.rds` | PNG images | Produces gene expression feature plots |
| **plot_cluster_composition.R** | R | Seurat `.rds` | Bar plots | Displays cell type composition by condition/sample |
| **interactive_visualization_app.R** | R | Seurat `.rds` | Shiny web app | Launches an interactive viewer for metadata and gene expression |
| **generate_visualization_report.R** | R | Directory of plots | PDF | Combines generated figures into a report |
| **.gitkeep** | text | — | — | Keeps folder versioned when empty |

---

## ⚙️ Usage Examples

### 1️⃣ Embedding Visualization
```bash
projects/transcriptomics/sc_rna/sc_visualization/plot_embeddings.R \
  -i results/annotated_seurat.rds \
  -m cell_type \
  -o results/visualization/embeddings
```

### 2️⃣ Feature Expression Plots
```bash
projects/transcriptomics/sc_rna/sc_visualization/plot_features.R \
  -i results/annotated_seurat.rds \
  -g CD3D,MS4A1,LYZ \
  -o results/visualization/features
```

### 3️⃣ Cluster Composition
```bash
projects/transcriptomics/sc_rna/sc_visualization/plot_cluster_composition.R \
  -i results/annotated_seurat.rds \
  -g cell_type \
  -c condition \
  -o results/visualization/composition
```

### 4️⃣ Interactive Exploration (Shiny)
```bash
projects/transcriptomics/sc_rna/sc_visualization/interactive_visualization_app.R \
  -i results/annotated_seurat.rds
```
Then open the provided URL (default: `http://127.0.0.1:xxxx`) in your browser.

### 5️⃣ Generate Summary Report
```bash
projects/transcriptomics/sc_rna/sc_visualization/generate_visualization_report.R \
  -i results/visualization/embeddings \
  -o reports/visualization_summary.pdf
```

---

## 🧠 Best Practices

- Use **normalized and annotated** Seurat/AnnData objects as input.  
- For UMAP/tSNE, ensure dimensional reductions are computed (`RunUMAP`, `RunTSNE`).  
- Always include clear legends and color schemes consistent across figures.  
- Limit the number of genes per `FeaturePlot` for readability.  
- For cluster composition, standardize metadata field names (`cell_type`, `condition`).  
- Save both PNG and PDF versions for publication and further editing.  
- For Shiny apps, run locally or deploy via **Shiny Server / RStudio Connect** for team exploration.  
- Maintain a consistent plotting theme across the entire project (e.g., `theme_minimal()` or `theme_classic()`).

---

## 🧾 Attribution

Created and maintained by **Farouk Saaidia (2025)**.  
For research or educational reuse, please cite:

> Saaidia F. (2025). *Single-Cell RNA-seq Visualization and Interactive Exploration Module.*

