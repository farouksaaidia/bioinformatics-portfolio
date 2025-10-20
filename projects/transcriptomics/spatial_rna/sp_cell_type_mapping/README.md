# sp_cell_type_mapping

Cell-type mapping and label transfer module for spatial transcriptomics.  
Integrates scRNA-seq reference information with spatial expression data to infer spatial distributions of cell types.

---

## 📂 Folder Purpose

`projects/transcriptomics/spatial_rna/sp_cell_type_mapping/`  
Implements multiple strategies for mapping reference cell-type annotations to spatial transcriptomics datasets, including anchor-based transfer, probabilistic mapping, and model comparison.

---

## 🧩 Scripts Overview

| Script | Language | Input | Output | Description |
|--------|-----------|--------|---------|-------------|
| **run_label_transfer_seurat.R** | R | scRNA `.rds`, spatial `.rds` | annotated `.rds` | Seurat-based label transfer |
| **run_tangram_mapping.py** | Python | scRNA `.h5ad`, spatial `.h5ad` | annotated `.h5ad` | Tangram probabilistic mapping |
| **run_cell2location_mapping.py** | Python | scRNA `.h5ad`, spatial `.h5ad` | annotated `.h5ad` | cell2location Bayesian model |
| **compare_mapping_methods.py** | Python | two CSVs | metrics + confusion matrix | Compare mappings across methods |
| **generate_mapping_report.R** | R | annotated `.rds` | PDF report | Visual overlay of predicted cell types |
| **.gitkeep** | text | — | — | Keeps folder tracked by Git |

---

## ⚙️ Example Workflows

### 1️⃣ Seurat label transfer
Rscript projects/transcriptomics/spatial_rna/sp_cell_type_mapping/run_label_transfer_seurat.R \
  -r results/sc_reference/reference_immune.rds \
  -q results/spatial/visium_sample.rds \
  -o results/mapping/visium_label_transfer.rds

### 2️⃣ Tangram probabilistic mapping
projects/transcriptomics/spatial_rna/sp_cell_type_mapping/run_tangram_mapping.py \
  -r results/sc_reference/reference.h5ad \
  -s results/spatial/sample01.h5ad \
  -o results/mapping/sample01_tangram.h5ad

### 3️⃣ cell2location Bayesian inference
projects/transcriptomics/spatial_rna/sp_cell_type_mapping/run_cell2location_mapping.py \
  -r results/sc_reference/reference.h5ad \
  -s results/spatial/sample01.h5ad \
  -o results/mapping/sample01_cell2location.h5ad

### 4️⃣ Compare mapping outputs
projects/transcriptomics/spatial_rna/sp_cell_type_mapping/compare_mapping_methods.py \
  --a results/mapping/seurat_pred.csv \
  --b results/mapping/tangram_pred.csv \
  --out_prefix results/mapping/comparison

### 5️⃣ Generate visual mapping report
Rscript projects/transcriptomics/spatial_rna/sp_cell_type_mapping/generate_mapping_report.R \
  -i results/mapping/visium_label_transfer.rds \
  -o reports/mapping_report.pdf

---

## 🧠 Best Practices

- Use high-quality **scRNA-seq reference** with accurate cell-type annotations.  
- Ensure gene overlap between reference and spatial datasets.  
- Tangram and cell2location perform best with well-normalized data.  
- Compare multiple mapping methods for consistency — no single approach is perfect.  
- Visual inspection of mappings is essential for spatial context validation.  

---

## 🧾 Attribution

Created and maintained by **Farouk Saaidia (2025)**.
