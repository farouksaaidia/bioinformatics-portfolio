# sp_preprocessing

Preprocessing utilities for spatial transcriptomics (Visium / Slide-seq / Stereo-seq).  
This module focuses on loading spatial outputs, integrating histology images, computing spot-level QC, and exporting interoperable Seurat/AnnData objects.

---

## 📂 Folder Purpose

`projects/transcriptomics/spatial_rna/sp_preprocessing/`  
Contains all tools required to load raw spatial transcriptomics data, attach histology images, compute spatial quality metrics, and prepare interoperable analysis objects (`.rds`, `.h5seurat`, `.h5ad`).

---

## 🧩 Scripts Overview

| Script | Language | Input | Output | Description |
|--------|-----------|--------|---------|-------------|
| **load_visium_seurat.R** | R | 10x Visium sample folder | Seurat `.rds` | Robust loader using `Load10X_Spatial`, adds sample metadata and QC |
| **load_visium_scanpy.py** | Python | 10x Visium sample folder | AnnData `.h5ad` | Loads spatial counts + image metadata with Scanpy `read_visium()` |
| **integrate_histology_image.R** | R | Seurat `.rds` + image (PNG/JPEG) | updated Seurat `.rds` | Integrates custom histology image and computes image-derived intensity metrics per spot |
| **extract_spatial_qc_metrics.R** | R | Seurat `.rds` | spot QC CSV + updated `.rds` | Computes counts, features, %MT, and spatial neighbor metrics |
| **export_standard_objects.sh** | Bash + R | Seurat `.rds` | `.h5seurat` + `.h5ad` | Exports canonical Seurat object to interoperable formats via SeuratDisk |

---

## ⚙️ Usage Examples

### 1️⃣ Load Visium (Seurat)
```bash
Rscript projects/transcriptomics/spatial_rna/sp_preprocessing/load_visium_seurat.R \
  --sample_dir data/sample1/ --sample_name SAMPLE01 \
  --output results/preprocessing/SAMPLE01_seurat.rds
```

### 2️⃣ Load Visium (Scanpy / AnnData)
```bash
projects/transcriptomics/spatial_rna/sp_preprocessing/load_visium_scanpy.py \
  -s data/sample1/ -n SAMPLE01 -o results/preprocessing/SAMPLE01.h5ad
```

### 3️⃣ Integrate Custom Histology Image
```bash
Rscript projects/transcriptomics/spatial_rna/sp_preprocessing/integrate_histology_image.R \
  -i results/preprocessing/SAMPLE01_seurat.rds \
  -img data/images/SAMPLE01_hires.png \
  -o results/preprocessing/SAMPLE01_with_image.rds
```

### 4️⃣ Extract Spatial QC Metrics
```bash
Rscript projects/transcriptomics/spatial_rna/sp_preprocessing/extract_spatial_qc_metrics.R \
  -i results/preprocessing/SAMPLE01_with_image.rds \
  -o results/qc/SAMPLE01_qc
```

### 5️⃣ Export for Python Workflows
```bash
projects/transcriptomics/spatial_rna/sp_preprocessing/export_standard_objects.sh \
  -r results/preprocessing/SAMPLE01_with_image.rds \
  -o results/interoperability/SAMPLE01
```

---

## 🧠 Best Practices

- Prefer **Space Ranger outputs** for compatibility with `Load10X_Spatial()`.  
- Store histology images in **lossless PNG** when possible to preserve color metrics.  
- Keep consistent metadata fields: `sample_id`, `percent_mt`, and `histology_mean_intensity`.  
- Always verify image alignment after loading — coordinate mismatches are a common pitfall.  
- Use `SeuratDisk` for robust Seurat ↔ AnnData interoperability.  
- For reproducibility, record R/Python package versions via `sessionInfo()` or `pip freeze`.

---

## 🧾 Attribution

Created and maintained by **Farouk Saaidia (2025)**.  
For research or educational reuse, please cite:

> Saaidia F. (2025). *Spatial Transcriptomics Preprocessing Module.*

