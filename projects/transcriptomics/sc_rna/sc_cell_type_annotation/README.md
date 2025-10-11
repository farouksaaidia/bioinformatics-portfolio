# sc_cell_type_annotation

Comprehensive single-cell RNA-seq annotation module.  
This directory contains all scripts required to perform cell-type annotation, comparison, validation, and reporting for single-cell RNA-seq data.  
It integrates **marker-based**, **automated (CellTypist/Seurat Transfer)**, **ensemble**, and **ontology-aware** annotation workflows.

---

## 📂 Folder Purpose

`projects/transcriptomics/sc_rna/sc_cell_type_annotation/`  
Centralized module for annotation workflows in single-cell transcriptomics pipelines.  
All scripts are self-contained, CLI-based, and production-ready.

---

## 🧩 Scripts Overview

| Script | Language | Input | Output | Description |
|--------|-----------|--------|---------|-------------|
| **annotate_clusters_marker_based.R** | R | clustered Seurat `.rds`, markers CSV | annotated Seurat `.rds` | Marker-based cluster annotation using known gene markers and enrichment logic |
| **annotate_clusters_automated_scanpy.py** | Python | clustered `.h5ad` | annotated `.h5ad` (+ optional CSV) | Automated annotation via CellTypist reference models |
| **seurat_transfer_label_transfer.R** | R | reference `.rds`, query `.rds` | annotated query `.rds` | High-quality label transfer between Seurat objects using `FindTransferAnchors` + `TransferData` |
| **ensemble_annotation.py** | Python | `.h5ad` + optional marker/SingleR CSVs | `.h5ad` + votes CSV | Combine multiple annotations (CellTypist, marker-based, SingleR) into a consensus label with confidence scoring |
| **ontology_map_hierarchy.R** | R | annotated Seurat `.rds` | `.rds` with ontology metadata | Maps cell-type labels to Cell Ontology IDs and hierarchical parent relationships |
| **marker_enrichment_stats.R** | R | Seurat `.rds`, marker sets CSV | CSVs + RDS in output dir | Marker enrichment analysis (Fisher exact + AUCell) per cluster and per cell |
| **benchmark_annotations.py** | Python | gold CSV + predicted CSV | metrics CSVs + confusion PNG | Benchmark predicted annotations using precision/recall/F1 and confusion matrices |
| **compare_annotation_methods.py** | Python | two annotation CSVs | confusion matrix + metrics | Quantitatively compare results between two annotation methods (e.g., marker vs automated) |
| **generate_annotation_report.R** | R | annotated Seurat `.rds` | PDF + CSV outputs | Summary plots (UMAP colored by annotation + cell-type counts) |
| **interactive_annotation_report.R** | R | annotated Seurat `.rds` | interactive HTML report | Interactive report with UMAP, confusion heatmap, and label summaries via RMarkdown |
| **export_annotated_objects.sh** | Bash | annotated `.rds` / `.h5ad` | exported files + manifest | Export and organize final annotated Seurat and Scanpy objects for downstream use |
| **.gitkeep** | text | — | — | Keeps directory structure under version control |

---

## ⚙️ Usage Examples

### 1️⃣ Marker-based Annotation
```bash
projects/transcriptomics/sc_rna/sc_cell_type_annotation/annotate_clusters_marker_based.R \
  -i data/clustered_seurat.rds \
  -m data/known_markers.csv \
  -o results/annotated_marker_based.rds
```

### 2️⃣ Automated Annotation (CellTypist)
```bash
projects/transcriptomics/sc_rna/sc_cell_type_annotation/annotate_clusters_automated_scanpy.py \
  -i data/clustered.h5ad \
  -m Immune_All_Low.pkl \
  -o results/annotated_celltypist.h5ad \
  --pred_csv results/predictions.csv
```

### 3️⃣ Seurat Label Transfer
```bash
projects/transcriptomics/sc_rna/sc_cell_type_annotation/seurat_transfer_label_transfer.R \
  -r reference/ref_seurat.rds \
  -q query/query_seurat.rds \
  -o results/annotated_query_transfer.rds
```

### 4️⃣ Ensemble Consensus Annotation
```bash
projects/transcriptomics/sc_rna/sc_cell_type_annotation/ensemble_annotation.py \
  -i results/annotated_celltypist.h5ad \
  --marker_csv results/marker_annotation.csv \
  --singler_csv results/singler_annotation.csv \
  -o results/ensemble_annotation.h5ad \
  --out_votes results/ensemble_votes.csv
```

### 5️⃣ Ontology Mapping
```bash
projects/transcriptomics/sc_rna/sc_cell_type_annotation/ontology_map_hierarchy.R \
  -i results/ensemble_annotation.rds \
  -o results/ensemble_annotation_ontology.rds
```

### 6️⃣ Marker Enrichment Statistics
```bash
projects/transcriptomics/sc_rna/sc_cell_type_annotation/marker_enrichment_stats.R \
  -i results/ensemble_annotation.rds \
  -m data/marker_sets.csv \
  -o results/marker_enrichment
```

### 7️⃣ Benchmarking Annotations
```bash
projects/transcriptomics/sc_rna/sc_cell_type_annotation/benchmark_annotations.py \
  --gold data/gold_labels.csv \
  --pred results/predicted_labels.csv \
  --out_prefix results/benchmark
```

### 8️⃣ Compare Annotation Methods
```bash
projects/transcriptomics/sc_rna/sc_cell_type_annotation/compare_annotation_methods.py \
  -a results/marker_based_annotations.csv \
  -b results/automated_annotations.csv \
  -o results/comparison
```

### 9️⃣ Generate Summary Report (Static PDF)
```bash
projects/transcriptomics/sc_rna/sc_cell_type_annotation/generate_annotation_report.R \
  -i results/ensemble_annotation.rds \
  -o reports/static_report
```

### 🔟 Interactive HTML Report
```bash
projects/transcriptomics/sc_rna/sc_cell_type_annotation/interactive_annotation_report.R \
  -i results/ensemble_annotation.rds \
  -o reports/interactive_report.html
```

### 1️⃣1️⃣ Export Annotated Objects
```bash
projects/transcriptomics/sc_rna/sc_cell_type_annotation/export_annotated_objects.sh \
  -r results/ensemble_annotation.rds \
  -h results/annotated_celltypist.h5ad \
  -o exports/final_outputs
```

---

## 🧠 Best Practices

- Always normalize and log-transform Seurat objects before label transfer or marker detection.  
- Maintain consistent feature spaces between reference and query datasets.  
- Use standardized metadata fields:
  - `seurat_clusters`
  - `predicted_cell_type`
  - `predicted_cell_type_transfer`
  - `ensemble_label`
- Use HGNC symbols for gene identifiers in all CSVs and reference tables.  
- Record session/package versions for reproducibility:
  - In R: `sessionInfo()`
  - In Python: `pip freeze` or `conda env export`
- Use Git LFS or cloud storage for large `.rds` or `.h5ad` files.  
- Extend ontology mappings (`ontology_map_hierarchy.R`) with full Cell Ontology references for production.  
- Prefer benchmark datasets with expert-curated cell-type annotations for evaluation.

---

## 🧾 Attribution
Created and maintained by **Farouk Saaidia (2025)**.  
For academic or collaborative use, please cite this module as:

> Saaidia F. (2025). *Single-cell RNA annotation toolkit – sc_cell_type_annotation module.*

