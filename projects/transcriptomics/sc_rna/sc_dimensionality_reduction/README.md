# 🧬 Single-Cell RNA-seq: Dimensionality Reduction

## 📖 Overview
After normalization, **dimensionality reduction** is one of the most important steps in single-cell RNA-seq (scRNA-seq) analysis.  
It transforms thousands of gene expression measurements per cell into a few informative components that capture biological variability while reducing noise.

This step is critical for:
- Identifying major sources of variation (e.g., cell type, cell cycle, technical batch)
- Improving visualization and clustering
- Enabling trajectory inference and differential expression downstream

---

## 📂 Included Scripts

| Script | Framework | Input | Output | When to Use | Particularity |
|---------|------------|--------|----------|---------------|----------------|
| **run_pca_seurat.R** | Seurat (R) | Normalized Seurat object (`.rds`) | Seurat object with PCA embeddings | Always — first linear DR step | Fast, interpretable; prepares data for UMAP/t-SNE |
| **run_pca_scanpy.py** | Scanpy (Python) | AnnData object (`.h5ad`) | AnnData with PCA embeddings | Always — before non-linear embedding | Standardized workflow in Scanpy; integrates well with batch correction |
| **run_umap_seurat.R** | Seurat (R) | PCA-reduced Seurat object | UMAP 2D coordinates and plots | For visualization and clustering | Non-linear embedding emphasizing global structure |
| **run_umap_scanpy.py** | Scanpy (Python) | PCA-reduced AnnData | UMAP coordinates and `.png` plots | For visualization and downstream annotation | Fast, robust; integrates with Scanpy pipeline |
| **run_tsne.R** | Seurat (R) | PCA-reduced Seurat object | t-SNE embeddings and `.png` plots | For local neighborhood visualization | Alternative to UMAP; highlights subtle cell-state variation |
| **visualize_embeddings.R** | Seurat (R) | Seurat object with embeddings | Multi-panel `.png` or `.pdf` plots | After UMAP/t-SNE | Combines QC metrics and gene expression on embedding plots |

---

## 🧠 Notes & Best Practices

- Always perform **PCA first**, even when using UMAP or t-SNE.  
  These non-linear methods rely on PCA-reduced features as input.
- Use **UMAP** for general visualization — it better preserves global topology than t-SNE.  
- **t-SNE** can be valuable to explore fine-grained subclusters or rare cell states.
- Dimensionality reduction also **removes technical noise** and facilitates clustering and trajectory inference.
- Check the **explained variance** of your PCA to ensure enough components are kept (usually 30–50 PCs).
- Visualize and compare results from Seurat and Scanpy to confirm consistency.

---

📦 **Output summary**
- `.rds` and `.h5ad` objects containing reduced dimensions
- `.png` and `.pdf` plots of UMAP/t-SNE embeddings
- Summary tables of explained variance
