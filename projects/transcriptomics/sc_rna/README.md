# 🧬 sc_rna — Single-Cell RNA-Seq Analysis Pipeline

Comprehensive, modular pipeline for **single-cell RNA sequencing (scRNA-seq)** analysis.  
Each subfolder corresponds to a major analysis stage — from raw data preprocessing to trajectory inference and visualization.  
The framework supports both **R (Seurat, Monocle)** and **Python (Scanpy)** backends, enabling hybrid and cross-validated workflows.

---

## 🧭 Overview of Pipeline Stages

| Stage | Folder | Goal | Key Tools |
|-------|---------|------|-----------|
| **1️⃣ Preprocessing** | `sc_preprocessing/` | Generate count matrices from raw FASTQ files using standard mappers | CellRanger, STARsolo, Alevin, DropSeqTools, Kallisto-Bustools, ZUMIs |
| **2️⃣ Quality Control** | `sc_qc/` | Filter poor-quality cells, remove doublets, assess sequencing depth & mitochondrial content | Seurat, Scanpy, Scrublet, EmptyDrops, FastQC, MultiQC |
| **3️⃣ Normalization & Batch Correction** | `sc_normalization/` | Normalize, scale, and correct batch effects | Seurat SCTransform, scran, log-normalization, MNN, Harmony |
| **4️⃣ Dimensionality Reduction** | `sc_dimensionality_reduction/` | Reduce dimensionality for visualization and clustering | PCA, UMAP, t-SNE, Scanpy embeddings |
| **5️⃣ Clustering** | `sc_clustering/` | Identify cell populations using resolution-optimized clustering | Seurat/Scanpy Leiden-Louvain clustering, resolution tuning |
| **6️⃣ Cell Type Annotation** | `sc_cell_type_annotation/` | Assign biological identities using markers and automated models | Marker-based, Seurat label transfer, CellTypist, SingleR, ensemble annotation |
| **7️⃣ Differential Expression** | `sc_differential_expression/` | Identify genes driving differences across clusters or conditions | Seurat FindMarkers, tradeSeq, GO/KEGG enrichment |
| **8️⃣ Trajectory Inference** | `sc_trajectory_inference/` | Infer pseudotime and lineage progression | Monocle3, Slingshot, Scanpy DPT, tradeSeq |
| **9️⃣ Visualization** | `sc_visualization/` | Generate publication-ready static and interactive plots | Seurat, ggplot2, Shiny, custom R/Python visualization scripts |

---

## 🧩 Folder-by-Folder Summary

### 🧱 `sc_preprocessing/`
Prepares gene expression matrices from raw reads using common single-cell pipelines.  
**Includes:**  
- `run_cellranger.sh`, `run_starsolo.sh`, `run_kallisto_bustools.sh`, `run_alevin.sh`, `run_dropseqtools.sh`, `run_zumis.sh`  
**Output:** filtered count matrices in `.h5`, `.mtx`, or `.loom`.

---

### 🔍 `sc_qc/`
Performs rigorous QC to retain high-quality cells and genes.  
**Includes:**  
- `run_preprocess_seurat.R`, `run_preprocess_scanpy.py`, `run_scrublet_scanpy.py`, `run_emptydrops.R`, `generate_qc_report.py`  
**Output:** QC plots, filtered objects, and metrics summary.

---

### ⚖️ `sc_normalization/`
Handles normalization, scaling, and batch correction across datasets.  
**Includes:**  
- `run_log_normalization.R`, `run_sctransform.R`, `run_scran_pooling.R`, `run_multi_sample_normalization.R`, `run_batch_correction.R`, `run_scaling.R`  
**Output:** normalized and integrated Seurat/Scanpy objects.

---

### 🔢 `sc_dimensionality_reduction/`
Computes PCA, UMAP, and t-SNE embeddings for downstream clustering and visualization.  
**Includes:**  
- `run_pca_seurat.R`, `run_umap_scanpy.py`, `run_tsne.R`, `visualize_embeddings.R`  
**Output:** objects with dimensional reductions stored and plots of embeddings.

---

### 🧬 `sc_clustering/`
Groups cells into transcriptionally similar populations.  
**Includes:**  
- `run_clustering_seurat.R`, `run_clustering_scanpy.py`, `optimize_resolution.R`, `compare_clusterings.py`, `export_cluster_markers.R`  
**Output:** cluster-assigned objects and marker gene tables.

---

### 🧫 `sc_cell_type_annotation/`
Assigns biological identities to clusters using multiple approaches.  
**Includes:**  
- `annotate_clusters_marker_based.R`, `annotate_clusters_automated_scanpy.py`, `seurat_transfer_label_transfer.R`, `ensemble_annotation.py`, `ontology_map_hierarchy.R`, `marker_enrichment_stats.R`  
- `generate_annotation_report.R`, `interactive_annotation_report.R`, `benchmark_annotations.py`, `compare_annotation_methods.py`, `export_annotated_objects.sh`  
**Output:** annotated objects, summary plots, and interactive HTML reports.

---

### 🧠 `sc_differential_expression/`
Identifies genes and pathways differentiating clusters or conditions.  
**Includes:**  
- `find_differential_genes.R`, `visualize_deg_results.R`, `pathway_enrichment_analysis.R`, `compare_deg_sets.py`, `generate_deg_report.R`  
**Output:** DEG tables, enrichment results, and PDF summaries.

---

### 🧩 `sc_trajectory_inference/`
Reconstructs lineage progression and pseudotime.  
**Includes:**  
- `run_monocle3.R`, `run_slingshot.R`, `run_scanpy_dpt.py`, `pseudotime_deg.R`, `visualize_trajectory.R`, `generate_trajectory_report.R`  
**Output:** pseudotime values, DEGs along trajectory, lineage visualizations.

---

### 🎨 `sc_visualization/`
Centralized visualization utilities and dashboards.  
**Includes:**  
- `plot_embeddings.R`, `plot_cluster_composition.R`, `plot_features.R`, `interactive_visualization_app.R`, `generate_visualization_report.R`  
**Output:** static and interactive visual summaries across all modules.

---

## 🔗 Full Pipeline Flow

```
Raw FASTQ
   │
   ├──▶ sc_preprocessing (CellRanger / STARsolo / etc.)
   │
   ├──▶ sc_qc (Filter, Scrublet, QC report)
   │
   ├──▶ sc_normalization (SCTransform / scran / Harmony)
   │
   ├──▶ sc_dimensionality_reduction (PCA → UMAP)
   │
   ├──▶ sc_clustering (Leiden/Louvain)
   │
   ├──▶ sc_cell_type_annotation (Marker / CellTypist / Seurat transfer)
   │
   ├──▶ sc_differential_expression (DEG, GO/KEGG)
   │
   ├──▶ sc_trajectory_inference (Monocle3 / Slingshot / DPT)
   │
   └──▶ sc_visualization (Final static + interactive plots)
```

---

## 🧠 Best Practices

- Keep consistent file naming: `{sample_id}_{stage}.rds` or `.h5ad`.  
- Always validate clustering before annotation (e.g., resolution optimization).  
- Integrate both R and Python backends for cross-validation of results.  
- Maintain reproducibility via environment files (`renv.lock`, `requirements.txt`).  
- Commit changes using semantic messages (e.g., `feat`, `fix`, `docs`, `chore`).  
- Run periodic `make clean` or `snakemake --cleanup-metadata` to reduce artifacts.

---

## 🧾 Attribution
Developed and maintained by **Farouk Saaidia (2025)**.  
For research reuse, please cite:

> Saaidia F. (2025). *Modular Single-Cell RNA-Seq Pipeline for Comprehensive Transcriptomic Analysis.*

