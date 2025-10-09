# Single-Cell RNA-seq — Clustering Module

This module performs **cell clustering** on normalized and dimensionally reduced single-cell RNA-seq data using Seurat (R) and Scanpy (Python). Clustering identifies transcriptionally distinct cell populations, which is crucial for uncovering cell types and heterogeneity in complex tissues.

## Scripts Overview

| Script | Language | Input | Output | When to Use | Key Features |
|--------|----------|-------|--------|-------------|--------------|
| run_clustering_seurat.R | R | Normalized `.rds` (Seurat object) | Clustered `.rds` with `seurat_clusters` metadata | Standard R-based clustering | Louvain/Leiden algorithms; single & batch sample support |
| run_clustering_scanpy.py | Python | Normalized `.h5ad` (Scanpy object) | Clustered `.h5ad` with cluster metadata | Standard Python-based clustering | Supports multi-resolution testing; parallelizable |
| optimize_resolution.R | R | Clustered `.rds` object | `resolution_cluster_counts.csv`, `resolution_vs_clusters.png` | Tune clustering resolution | Automated resolution sweep; visual summary |
| compare_clusterings.py | Python | Two `.h5ad` or `.rds` cluster-labeled files | Contingency table CSV, ARI score, heatmap | Compare methods or parameter choices | Computes Adjusted Rand Index and visual comparison |
| export_cluster_markers.R | R | Clustered `.rds` object | `markers_all_clusters.csv`, per-cluster `top_markers.csv` | For downstream annotation | Finds marker genes using `FindAllMarkers()` in Seurat |

## Notes
- Graph construction: `FindNeighbors()` (Seurat) / `sc.pp.neighbors()` (Scanpy)  
- Community detection: Louvain and Leiden algorithms  
- Resolution parameter: Controls granularity; higher = finer subclusters  
- Automated resolution sweep: Helps identify biologically meaningful clusters  
- Marker export: Facilitates downstream cell-type annotation and functional interpretation  

