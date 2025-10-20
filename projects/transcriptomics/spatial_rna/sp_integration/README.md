# sp_integration

Integration module for spatial transcriptomics: merging slides, batch correction, and joint integration with single-cell references.  
This module supports merging Seurat spatial objects, Harmony-based integration, Scanpy batch-correction (BBKNN/Scanorama), anchor-based Seurat label transfer, and QC reporting.

---

## 📂 Folder Purpose

`projects/transcriptomics/spatial_rna/sp_integration/`  
Prepare integrated multi-sample spatial datasets for joint downstream analysis (clustering, annotation, DE, inference) and transfer labels from high-quality single-cell references.

---

## 🧩 Scripts Overview

| Script | Language | Input | Output | Description |
|--------|-----------|--------|---------|-------------|
| **merge_spatial_datasets.R** | R | multiple Seurat `.rds` | merged Seurat `.rds` | Merge slides preserving images and metadata; prefixes cell/spot names |
| **integrate_spatial_seurat_harmony.R** | R | merged or multiple Seurat `.rds` | integrated Seurat `.rds` | Run SCTransform/PCA then Harmony on PCs to remove batch/sample effects |
| **batch_correction_scanpy.py** | Python | AnnData `.h5ad` | integrated `.h5ad` | BBKNN or Scanorama options for Python-based batch correction |
| **integrate_spatial_with_scrna.R** | R | scRNA ref `.rds`, spatial `.rds` | annotated spatial `.rds` | Seurat anchor-based integration and label transfer from scRNA reference |
| **generate_integration_report.R** | R | integrated Seurat `.rds` | PDF report | Embedding visualizations, batch mixing checks, and label-transfer QC |
| **.gitkeep** | text | — | — | Keeps folder tracked by Git |

---

## ⚙️ Example Workflows

### 1️⃣ Merge multiple spatial Seurat objects
Rscript projects/transcriptomics/spatial_rna/sp_integration/merge_spatial_datasets.R \
  -i "sample1.rds,sample2.rds" -n "S1,S2" -o results/integration/merged_samples.rds

### 2️⃣ Integrate merged slides with Harmony
Rscript projects/transcriptomics/spatial_rna/sp_integration/integrate_spatial_seurat_harmony.R \
  -i results/integration/merged_samples.rds -b sample_id -o results/integration/merged_integrated.rds

### 3️⃣ Batch-correct using Scanpy (BBKNN)
projects/transcriptomics/spatial_rna/sp_integration/batch_correction_scanpy.py \
  -i results/integration/merged_samples.h5ad -o results/integration/merged_bbknn.h5ad \
  --batch_key sample_id --method bbknn

### 4️⃣ Integrate spatial with scRNA reference and transfer labels
Rscript projects/transcriptomics/spatial_rna/sp_integration/integrate_spatial_with_scrna.R \
  -r references/sc_reference.rds -s results/preprocessing/SAMPLE01_seurat.rds \
  -o results/integration/SAMPLE01_label_transfer.rds

### 5️⃣ Generate integration QC report
Rscript projects/transcriptomics/spatial_rna/sp_integration/generate_integration_report.R \
  -i results/integration/merged_integrated.rds -b sample_id -o reports/integration_report.pdf

---

## 🧠 Notes & Best Practices

- Prefer SCTransform for Seurat-based integration; keep normalized assays for comparisons.  
- Harmony operates on PCs — make sure PCs capture biological variation before running Harmony.  
- For Python workflows, BBKNN is fast and robust for batch mixing; Scanorama is an alternative for more aggressive correction.  
- Always inspect cluster composition by batch/sample to ensure biological signals are preserved.  
- When transferring labels from scRNA references, verify gene overlap and normalization strategies.  
- Store both integrated embeddings and raw per-sample embeddings for reproducibility and diagnostics.

---

## 🧾 Attribution

Created and maintained by **Farouk Saaidia (2025)**.
