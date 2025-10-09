# 🧬 Single-Cell RNA-seq — Clustering Module

This module performs **cell clustering** on normalized and dimensionally reduced single-cell RNA-seq data.  
It uses **graph-based community detection** (Louvain, Leiden) through both **Seurat (R)** and **Scanpy (Python)** frameworks to identify transcriptionally distinct cell populations.  
Clustering is a critical step to uncover biological cell types and heterogeneity within complex tissues.

---

## 📁 Scripts Overview

| Script | Language | Input | Output | When to Use | Key Features |
|:--------|:----------:|:--------|:---------|:---------------|:----------------|
| `run_clustering_seurat.R` | R | Normalized `.rds` file (Seurat object) | Clustered `.rds` file with `seurat_clusters` metadata | Standard clustering in R pipelines | Implements Louvain/Leiden; supports batch or single-sample mode |
| `run_clustering_scanpy.py` | Python | Normalized `.h5ad` file (Scanpy object) | Clustered `.h5ad` with cluster metadata | Python-based workflows | Parallelizable; supports multi-resolution testing |
| `optimize_resolution.R` | R | Clustered `.rds` object | `resolution_cluster_counts.csv`, `resolution_vs_clusters.png` | To tune resolution for best biological separation | Automated resolution sweep with visual summary |
| `compare_clusterings.py` | Python | Two `.h5ad` or `.rds` files with cluster labels | Contingency table CSV, ARI score, heatmap | When comparing methods or parameters | Computes Adjusted Rand Index & visualization |
| `export_cluster_markers.R` | R | Clustered `.rds` object | `markers_all_clusters.csv`, per-cluster `top_markers.csv` | For downstream annotation & validation | Finds marker genes via Seurat’s `FindAllMarkers()` |

---

## ⚙️ Typical Workflow

1. **Input:** Normalized, PCA-reduced data (`.rds` or `.h5ad`)
2. **Run Clustering**
   ```bash
   # Seurat
   Rscript run_clustering_seurat.R -i input.rds -o clustered.rds -m leiden -r 0.6

   # Scanpy
   python run_clustering_scanpy.py -i input.h5ad -o clustered.h5ad -m louvain -r 0.5
Optimize Resolution

bash
Copy code
Rscript optimize_resolution.R -i clustered.rds -o results/ -r "0.2,0.4,0.6,0.8,1.0"
Compare Cluster Assignments

bash
Copy code
python compare_clusterings.py -i clustered.h5ad -c1 leiden_r0.4 -c2 louvain_r0.6 -o comparison/out
Export Marker Genes

bash
Copy code
Rscript export_cluster_markers.R -i clustered.rds -o markers/ -c seurat_clusters -p 0.05
🧩 Implementation Notes
Graph construction: FindNeighbors() (Seurat) / sc.pp.neighbors() (Scanpy)

Community detection: Louvain and Leiden algorithms

Resolution parameter: Controls granularity — higher values = finer subclustering

Automated optimization: Helps find biologically meaningful resolutions

Marker export: Enables interpretation and annotation downstream

🧠 Biological Insights
Discovers transcriptionally distinct cell populations within tissues

Reveals rare subpopulations and novel functional states

Forms the foundation for cell-type annotation, trajectory inference, and functional analysis

