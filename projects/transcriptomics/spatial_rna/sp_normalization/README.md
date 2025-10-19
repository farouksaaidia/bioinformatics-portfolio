# sp_normalization

Normalization, spatial smoothing, and cell-type deconvolution for spatial transcriptomics (Visium, Slide-seq, Stereo-seq).  
This module ensures consistent scaling, correction, and biological interpretability of spot-level gene expression before clustering or annotation.

---

## 📂 Folder Purpose

`projects/transcriptomics/spatial_rna/sp_normalization/`  
Contains scripts to normalize raw spatial counts (SCTransform or scran), enhance spatial resolution (BayesSpace), and deconvolve mixed cell populations (SPOTlight).

---

## 🧩 Scripts Overview

| Script | Language | Input | Output | Description |
|--------|-----------|--------|---------|-------------|
| **run_sctransform.R** | R | Seurat `.rds` | normalized Seurat `.rds` | SCTransform normalization with optional regression (percent_mt, batch, etc.) |
| **run_scran_pooling.R** | R | Seurat `.rds` | normalized Seurat `.rds` | scran pooled normalization (computeSumFactors, low UMI support) |
| **run_bayespace.R** | R | Seurat `.rds` (Visium) | updated Seurat `.rds` | BayesSpace spatial smoothing, MCMC-based enhancement, domain clustering |
| **run_spotlight_deconvolution.R** | R | scRNA reference `.rds` + spatial `.rds` | updated `.rds` | SPOTlight deconvolution to estimate cell-type proportions per spot |
| **.gitkeep** | text | — | — | Keeps folder tracked by Git |

---

## ⚙️ Example Usage

### 1️⃣ SCTransform normalization
Rscript projects/transcriptomics/spatial_rna/sp_normalization/run_sctransform.R \
  -i results/preprocessing/SAMPLE01_seurat.rds \
  -o results/normalized/SAMPLE01_sct.rds \
  --vars_to_regress percent_mt,batch

### 2️⃣ scran pooled normalization
Rscript projects/transcriptomics/spatial_rna/sp_normalization/run_scran_pooling.R \
  -i results/preprocessing/SAMPLE01_seurat.rds \
  -o results/normalized/SAMPLE01_scran.rds

### 3️⃣ BayesSpace spatial smoothing
Rscript projects/transcriptomics/spatial_rna/sp_normalization/run_bayespace.R \
  -i results/normalized/SAMPLE01_sct.rds \
  -k 7 -r 1000 \
  -o results/normalized/SAMPLE01_bayespace.rds

### 4️⃣ SPOTlight deconvolution
Rscript projects/transcriptomics/spatial_rna/sp_normalization/run_spotlight_deconvolution.R \
  -r references/sc_reference_seurat.rds \
  -s results/normalized/SAMPLE01_sct.rds \
  -l cell_type \
  -o results/normalized/SAMPLE01_spotlight.rds

---

## 🧠 Best Practices

- **SCTransform** is the default choice for Visium; regress technical variables like `percent_mt` or `batch`.  
- Use **scran** for lower-depth data (Slide-seq, Stereo-seq).  
- Run **BayesSpace** only on Visium with spatial coordinates; check domain maps visually.  
- Ensure single-cell reference for **SPOTlight** is well-annotated and harmonized.  
- Keep separate assays for each normalization (`SCT`, `scran`, `bayespace_enhanced`, `spotlight_props`).  
- Log software versions (Seurat, scran, SPOTlight, BayesSpace) for reproducibility.

---

## 🧾 Attribution

Created and maintained by **Farouk Saaidia (2025)**.
