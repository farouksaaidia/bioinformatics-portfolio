# sp_clustering

Spatial clustering and domain detection module for spatial transcriptomics datasets (Visium, Slide-seq, Stereo-seq).  
This stage identifies spatial domains using both graph-based clustering and Bayesian spatial modeling.

---

## 📂 Folder Purpose

`projects/transcriptomics/spatial_rna/sp_clustering/`  
Implements spatial domain detection, clustering comparison, and visualization across Seurat, Scanpy, and BayesSpace pipelines.

---

## 🧩 Scripts Overview

| Script | Language | Input | Output | Description |
|--------|-----------|--------|---------|-------------|
| **run_spatial_clustering_seurat.R** | R | Seurat `.rds` | clustered `.rds` | Graph-based clustering with optional spatial coordinate integration |
| **run_spatial_clustering_scanpy.py** | Python | `.h5ad` | clustered `.h5ad` | Leiden/Louvain clustering on spatial neighbor graph |
| **run_bayespace_domains.R** | R | Seurat `.rds` | `.rds` with BayesSpace domains | Bayesian domain detection via MCMC |
| **compare_spatial_clusterings.py** | Python | two CSVs | confusion matrix + ARI metrics | Compare two clusterings quantitatively |
| **generate_clustering_report.R** | R | Seurat `.rds` | PDF | Visual clustering report overlayed on tissue image |
| **.gitkeep** | text | — | — | Keeps folder tracked by Git |

---

## ⚙️ Example Workflows

### 1️⃣ Seurat spatial clustering
Rscript projects/transcriptomics/spatial_rna/sp_clustering/run_spatial_clustering_seurat.R \
  -i results/normalized/SAMPLE01_sct.rds \
  -o results/clustering/SAMPLE01_clusters.rds \
  --resolution 0.6 --use_spatial

### 2️⃣ Scanpy spatial clustering
projects/transcriptomics/spatial_rna/sp_clustering/run_spatial_clustering_scanpy.py \
  -i results/normalized/SAMPLE01.h5ad \
  -o results/clustering/SAMPLE01_leiden.h5ad \
  --method leiden --resolution 0.6

### 3️⃣ BayesSpace domain detection
Rscript projects/transcriptomics/spatial_rna/sp_clustering/run_bayespace_domains.R \
  -i results/normalized/SAMPLE01_bayespace.rds \
  -k 7 -r 1000 -o results/clustering/SAMPLE01_domains.rds

### 4️⃣ Compare clusterings
projects/transcriptomics/spatial_rna/sp_clustering/compare_spatial_clusterings.py \
  --a results/clustering/SAMPLE01_seurat_clusters.csv \
  --b results/clustering/SAMPLE01_bayespace_clusters.csv \
  --out_prefix results/clustering/comparison

### 5️⃣ Generate clustering report
Rscript projects/transcriptomics/spatial_rna/sp_clustering/generate_clustering_report.R \
  -i results/clustering/SAMPLE01_clusters.rds \
  -o reports/SAMPLE01_clustering_report.pdf

---

## 🧠 Best Practices

- **Match PCA/embedding basis** used for clustering across tools (same dimensional reductions).  
- Adjust resolution carefully: low values for broad tissue domains, high for fine subregions.  
- **BayesSpace** is most powerful for Visium data with histology alignment.  
- Always inspect cluster boundaries relative to tissue image — false boundaries may indicate low-quality regions.  
- Compare methods (Seurat vs BayesSpace vs Scanpy) using Adjusted Rand Index and spatial overlays.  

---

## 🧾 Attribution

Created and maintained by **Farouk Saaidia (2025)**.
