# sp_qc

Spatial Quality Control (QC) module for spatial transcriptomics datasets (Visium, Slide-seq, Stereo-seq).  
This step ensures high-quality data by computing per-spot metrics, detecting doublets, visualizing QC maps, and generating comprehensive PDF reports.

---

## 📂 Folder Purpose

`projects/transcriptomics/spatial_rna/sp_qc/`  
Contains all tools and scripts required to perform spatial QC, detect outliers or doublets, visualize QC indicators, and compile reports for downstream filtering and normalization.

---

## 🧩 Scripts Overview

| Script | Language | Input | Output | Description |
|--------|-----------|--------|---------|-------------|
| **run_spatial_qc_metrics.R** | R | Seurat `.rds` | `_qc_metrics.csv` + updated `.rds` | Compute counts, detected features, %MT, and spatial neighbor metrics |
| **detect_spatial_doublets_scanpy.py** | Python | AnnData `.h5ad` | `.h5ad` | Identify potential doublets using Scrublet and local density adjustment |
| **plot_spatial_qc_visuals.R** | R | Seurat `.rds` | PDF plots | Generate spatial maps for QC features (counts, features, mt%) |
| **generate_spatial_qc_report.R** | R | Seurat `.rds` + QC CSV | PDF report | Combine QC metrics and spatial plots into a summarized PDF report |

---

## ⚙️ Example Workflow

### 1️⃣ Compute QC Metrics
```bash
Rscript projects/transcriptomics/spatial_rna/sp_qc/run_spatial_qc_metrics.R \
  -i results/preprocessing/SAMPLE01_with_image.rds \
  -o results/qc/SAMPLE01
```

### 2️⃣ Detect Doublets
```bash
projects/transcriptomics/spatial_rna/sp_qc/detect_spatial_doublets_scanpy.py \
  -i results/preprocessing/SAMPLE01.h5ad \
  -o results/qc/SAMPLE01_doublets.h5ad
```

### 3️⃣ Generate QC Plots
```bash
Rscript projects/transcriptomics/spatial_rna/sp_qc/plot_spatial_qc_visuals.R \
  -i results/preprocessing/SAMPLE01_with_image.rds \
  -o results/qc/plots/
```

### 4️⃣ Create a Comprehensive QC Report
```bash
Rscript projects/transcriptomics/spatial_rna/sp_qc/generate_spatial_qc_report.R \
  -i results/preprocessing/SAMPLE01_with_image.rds \
  -q results/qc/SAMPLE01_qc_metrics.csv \
  -o reports/SAMPLE01_qc_report.pdf
```

---

## 🧠 Best Practices

- **Filter** low-quality or high-%MT spots before normalization.  
- **Cross-check** Scrublet doublets with image-based clusters or spot densities.  
- Use **UMAP/SpatialFeaturePlot** overlays to identify localized QC artifacts.  
- Maintain **consistent naming** for all intermediate outputs (`SAMPLEID_stage.qc_metrics.csv`).  
- Always **inspect histology alignment** before trusting spatial metrics.  
- Document **software versions** (`Seurat`, `Scanpy`, `Scrublet`, etc.) for reproducibility.

---

## 🧾 Attribution

Created and maintained by **Farouk Saaidia (2025)**.  
For research or educational reuse, please cite:

> Saaidia F. (2025). *Spatial Transcriptomics Quality Control and Doublet Detection Module.*

