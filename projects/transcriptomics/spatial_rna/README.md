# 🧬 Spatial Transcriptomics Analysis Pipeline (spatial_rna)

A modular, reproducible, and production-ready pipeline for analyzing spatial transcriptomics data — from raw Visium or Slide-seq files to annotated, integrated, and visualized tissue maps.  
Built using Seurat, Scanpy, BayesSpace, Giotto, Tangram, Cell2location, and Squidpy.

---

## 🚀 Overview

This workflow processes raw spatial transcriptomics data into interpretable maps of cell type organization and communication.  
It mirrors the `sc_rna` pipeline but adds spatial coordinates, histology alignment, and spatial inference capabilities.

---

## 🧭 Pipeline Flow

Raw spatial transcriptomics data (e.g., 10X Visium)

├──▶ sp_preprocessing → Load, format, and align raw data & histology

├──▶ sp_qc → Perform quality control on spots & coverage

├──▶ sp_normalization → Normalize counts, scale data, and deconvolve spots

├──▶ sp_dimensionality_reduction → Compute spatially aware PCA/UMAP embeddings

├──▶ sp_clustering → Identify tissue domains or regions

├──▶ sp_cell_type_mapping → Map scRNA-defined cell types onto spatial data

├──▶ sp_differential_expression → Identify region-specific markers & pathways

├──▶ sp_spatial_inference → Infer cell-cell communication and spatial signaling

├──▶ sp_integration → Integrate multiple slides or modalities (spatial + scRNA)

└──▶ sp_visualization → Generate final plots, interactive dashboards, and reports


---

## 🧩 Step-by-Step Modules

### 1️⃣ sp_preprocessing — Data Loading & Alignment
| Script | Description |
|--------|-------------|
| load_visium_seurat.R | Load 10x Visium data into a Seurat object with histology integration. |
| load_visium_scanpy.py | Load Visium data into an AnnData object for Scanpy workflows. |
| integrate_histology_image.R | Link histology images with spatial coordinates. |
| extract_spatial_qc_metrics.R | Compute nCount, nFeature, mitochondrial %, and spatial coverage. |
| export_standard_objects.sh | Export standardized `.rds` and `.h5ad` objects for downstream use. |

---

### 2️⃣ sp_qc — Spatial Quality Control
| Script | Description |
|--------|-------------|
| run_spatial_qc_metrics.R | Compute QC metrics per spot (nCount, nFeature, mito%). |
| detect_spatial_doublets_scanpy.py | Identify potential doublets spatially using Scanpy. |
| plot_spatial_qc_visuals.R | Plot QC maps over histology images. |
| generate_spatial_qc_report.R | Generate a full QC summary report (PDF). |

---

### 3️⃣ sp_normalization — Normalization & Deconvolution
| Script | Description |
|--------|-------------|
| run_sctransform.R | Apply variance-stabilizing normalization (SCTransform). |
| run_scran_pooling.R | Normalize using scran pooling for per-spot correction. |
| run_spotlight_deconvolution.R | Estimate cell-type proportions using scRNA reference. |
| run_bayespace.R | Apply BayesSpace normalization and spatial resolution enhancement. |

---

### 4️⃣ sp_dimensionality_reduction — Embedding Computation
| Script | Description |
|--------|-------------|
| run_spatial_pca.R | Perform spatially weighted PCA. |
| run_bayespace_pca.R | Compute spatially enhanced PCA using BayesSpace. |
| run_spatial_umap_scanpy.py | Generate UMAP embeddings considering spatial proximity. |
| generate_dimred_report.R | Summarize embeddings and create visualization reports. |

---

### 5️⃣ sp_clustering — Spatial Domain Detection
| Script | Description |
|--------|-------------|
| run_spatial_clustering_seurat.R | Cluster spatial spots in Seurat using expression and position. |
| run_spatial_clustering_scanpy.py | Cluster spatial data with Scanpy’s Leiden/Louvain. |
| run_bayespace_domains.R | Infer spatial domains using BayesSpace spatial priors. |
| compare_spatial_clusterings.py | Compare clustering outputs across algorithms. |
| generate_clustering_report.R | Generate domain summary plots and statistics. |

---

### 6️⃣ sp_cell_type_mapping — Cell-Type Mapping
| Script | Description |
|--------|-------------|
| run_label_transfer_seurat.R | Transfer scRNA-derived cell labels using Seurat anchors. |
| run_tangram_mapping.py | Map scRNA expression profiles to spatial coordinates using Tangram. |
| run_cell2location_mapping.py | Bayesian deconvolution of cell types with Cell2location. |
| compare_mapping_methods.py | Benchmark mapping agreement between methods. |
| generate_mapping_report.R | Produce visual summaries and confidence tables. |

---

### 7️⃣ sp_differential_expression — Region-Specific Markers
| Script | Description |
|--------|-------------|
| find_spatial_markers.R | Identify genes enriched in specific spatial domains. |
| differential_between_domains.R | Compare gene expression between selected regions. |
| spatial_moran_test.R | Compute Moran’s I for spatial autocorrelation. |
| pathway_enrichment_spatial.R | Run GO / KEGG enrichment on spatial DEGs. |
| generate_spatial_deg_report.R | Compile a comprehensive DE report with plots. |

---

### 8️⃣ sp_spatial_inference — Cell-Cell Communication
| Script | Description |
|--------|-------------|
| run_nichenet_inference.R | Predict ligand–receptor communication between domains. |
| run_giotto_neighborhood.R | Infer neighborhood-level cell proximity using Giotto. |
| run_squidpy_interactions.py | Compute spatial graph-based interactions via Squidpy. |
| run_liana_analysis.py | Integrate multi-method ligand–receptor predictions (LIANA). |
| compare_inference_methods.py | Evaluate concordance between inference methods. |
| generate_inference_report.R | Create visual and tabular summaries of inferred networks. |

---

### 9️⃣ sp_integration — Dataset Integration & Batch Correction
| Script | Description |
|--------|-------------|
| merge_spatial_datasets.R | Merge multiple spatial slides into one Seurat object. |
| integrate_spatial_seurat_harmony.R | Perform Harmony batch correction across slides. |
| integrate_spatial_with_scrna.R | Integrate scRNA and spatial datasets for multimodal mapping. |
| batch_correction_scanpy.py | Apply Scanpy-based batch correction. |
| generate_integration_report.R | Generate figures summarizing integration quality. |

---

### 🔟 sp_visualization — Visualization & Reporting
| Script | Description |
|--------|-------------|
| plot_spatial_features.R | Plot gene expression overlays on histology images. |
| plot_spatial_domains.R | Display spatial domains or clusters. |
| interactive_spatial_viewer.R | Launch an interactive Shiny app for spatial exploration. |
| generate_spatial_visual_report.R | Generate PDF reports combining spatial plots. |
| export_visual_assets.sh | Export figures, embeddings, and coordinates for publication. |

---

## 🧠 Best Practices
- Each step is standalone and callable via Snakemake or shell scripts.  
- Always validate QC metrics before normalization and clustering.  
- Use the same gene reference build across modules for reproducibility.  
- Typical multimodal order:  
  `sp_preprocessing → sp_qc → sp_normalization → sp_integration → sp_cell_type_mapping`.  
- Outputs from visualization modules are publication-ready and exportable as high-resolution figures.

---

## 🧾 Attribution
Developed by **Farouk Saaidia (2025)**.  
All scripts follow modular, version-controlled, and reproducible design for advanced spatial transcriptomics analysis.
