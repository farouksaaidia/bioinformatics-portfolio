# sp_dimensionality_reduction

Spatial-aware dimensionality reduction utilities for spatial transcriptomics (Visium, Slide-seq, Stereo-seq).  
This module computes PCA/UMAP embeddings that incorporate expression and — optionally — spatial context, and produces visual reports overlayed on tissue.

---

## 📂 Folder Purpose

`projects/transcriptomics/spatial_rna/sp_dimensionality_reduction/`  
Compute and store low-dimensional embeddings for spatial spots, including spatially-informed PCA and BayesSpace-enhanced PCA, and generate summary reports.

---

## 🧩 Scripts Overview

| Script | Language | Input | Output | Description |
|--------|-----------|--------|---------|-------------|
| **run_spatial_pca.R** | R | Seurat `.rds` | updated `.rds` with `pca` and optional `spatialPCA` reductions | Standard PCA on HVGs; optional simple spatial PCA by appending scaled XY coords as features |
| **run_spatial_umap_scanpy.py** | Python | AnnData `.h5ad` | `.h5ad` with PCA & UMAP | Scanpy PCA/neighbors/UMAP runner preserving `obsm['spatial']` |
| **run_bayespace_pca.R** | R | Seurat `.rds` (BayesSpace enhanced assay) | updated `.rds` with `bayesPCA` reduction | PCA on BayesSpace-enhanced assay to highlight spatial signals |
| **generate_dimred_report.R** | R | Seurat `.rds` | PDF report | Combined embedding plots and spatial overlays for QC and interpretation |
| **.gitkeep** | text | — | — | Keeps folder tracked by Git |

---

## ⚙️ Example Workflows

### 1️⃣ Standard PCA with optional spatial PCA
Rscript projects/transcriptomics/spatial_rna/sp_dimensionality_reduction/run_spatial_pca.R \
  -i results/normalized/SAMPLE01_sct.rds \
  -o results/dimred/SAMPLE01_pca.rds \
  --npcs 30 --spatial_pca

### 2️⃣ Scanpy PCA + UMAP (AnnData)
projects/transcriptomics/spatial_rna/sp_dimensionality_reduction/run_spatial_umap_scanpy.py \
  -i results/preprocessing/SAMPLE01.h5ad \
  -o results/dimred/SAMPLE01_umap.h5ad \
  --n_pcs 30 --n_neighbors 15

### 3️⃣ BayesSpace PCA (on enhanced assay)
Rscript projects/transcriptomics/spatial_rna/sp_dimensionality_reduction/run_bayespace_pca.R \
  -i results/normalized/SAMPLE01_bayespace.rds \
  -o results/dimred/SAMPLE01_bayespace_pca.rds

### 4️⃣ Generate combined embedding + spatial report
Rscript projects/transcriptomics/spatial_rna/sp_dimensionality_reduction/generate_dimred_report.R \
  -i results/dimred/SAMPLE01_pca.rds \
  -o reports/SAMPLE01_dimred_report.pdf \
  -f nFeature_RNA,percent_mt

---

## 🧠 Notes & Best Practices

- Spatial PCA in `run_spatial_pca.R` uses a pragmatic approach (appending scaled X/Y coordinates) — it is simple, reproducible, and frequently useful as a quick spatial embedding.  
- Use **BayesSpace-enhanced** assay PCA to highlight sub-spot structure produced by BayesSpace.  
- When using Scanpy for UMAP, ensure neighbors are computed with spatially relevant PCs if you want spatial preservation.  
- Visualize embeddings both as scatter plots and as overlays on histology to verify tissue-structure preservation.  
- Keep multiple reductions (`pca`, `spatialPCA`, `bayesPCA`, `umap`) simultaneously to compare results.

---

## 🧾 Attribution

Created and maintained by **Farouk Saaidia (2025)**.
