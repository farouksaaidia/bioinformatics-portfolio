# sc_trajectory_inference

Trajectory and pseudotime analysis module for single-cell RNA-seq.  
This module provides multiple complementary trajectory inference methods, pseudotime-based differential expression, visualization tools, and reporting.

---

## 📂 Folder Purpose

`projects/transcriptomics/sc_rna/sc_trajectory_inference/`  
Perform lineage inference, compute pseudotime, detect genes varying along trajectories, and produce visual summaries for interpretation.

---

## 🧩 Scripts Overview

| Script | Language | Input | Output | Description |
|--------|-----------|--------|---------|-------------|
| **run_monocle3.R** | R | Seurat `.rds` (normalized, PCA/UMAP) | Seurat `.rds` with `monocle3_pseudotime` | Monocle 3 graph learning and pseudotime ordering |
| **run_slingshot.R** | R | Seurat `.rds` | Seurat `.rds` with `sling_pseudotime` & slingshot object in `@misc` | Slingshot lineage inference on embedding |
| **run_scanpy_dpt.py** | Python | `.h5ad` | `.h5ad` with `dpt_pseudotime` | Diffusion pseudotime (Scanpy DPT) |
| **pseudotime_deg.R** | R | Seurat `.rds` with pseudotime | CSVs of genes changing along pseudotime | Uses tradeSeq for pseudotime-associated DE (GAMs) |
| **visualize_trajectory.R** | R | Seurat `.rds` | PNG plots (UMAP + pseudotime, gene trends) | Plot pseudotime on embedding and gene expression trends |
| **generate_trajectory_report.R** | R | Seurat `.rds` + optional DEG dir | PDF report | Combined summary of pseudotime analyses and top genes |
| **.gitkeep** | text | — | — | Keeps folder tracked under version control |

---

## ⚙️ Usage Examples

### 1️⃣ Monocle 3
```bash
projects/transcriptomics/sc_rna/sc_trajectory_inference/run_monocle3.R \
  -i results/annotated_seurat.rds \
  -o results/annotated_monocle3.rds
```

### 2️⃣ Slingshot
```bash
projects/transcriptomics/sc_rna/sc_trajectory_inference/run_slingshot.R \
  -i results/annotated_seurat.rds \
  -e umap \
  -o results/annotated_slingshot.rds
```

### 3️⃣ Scanpy Diffusion Pseudotime (DPT)
```bash
projects/transcriptomics/sc_rna/sc_trajectory_inference/run_scanpy_dpt.py \
  -i results/annotated.h5ad \
  -o results/annotated_dpt.h5ad
```

### 4️⃣ Pseudotime Differential Expression (tradeSeq)
```bash
projects/transcriptomics/sc_rna/sc_trajectory_inference/pseudotime_deg.R \
  -i results/annotated_monocle3.rds \
  -p monocle3_pseudotime \
  -o results/pseudotime_de
```

### 5️⃣ Visualization of Trajectories & Gene Trends
```bash
projects/transcriptomics/sc_rna/sc_trajectory_inference/visualize_trajectory.R \
  -i results/annotated_monocle3.rds \
  -p monocle3_pseudotime \
  -g GENE1,GENE2 \
  -o results/trajectory_plots
```

### 6️⃣ Summary Report (PDF)
```bash
projects/transcriptomics/sc_rna/sc_trajectory_inference/generate_trajectory_report.R \
  -i results/annotated_monocle3.rds \
  -d results/pseudotime_de \
  -o reports/trajectory_report.pdf
```

---

## 🧠 Best Practices

- Ensure Seurat or AnnData objects are normalized and contain dimensionality reductions (UMAP/PCA) before running trajectory inference.  
- Monocle3 works best with precomputed clustering and dimensionality reduction.  
- For branching trajectories, use tradeSeq for branch-aware DE analysis.  
- When using Scanpy DPT, ensure neighbors are computed (`sc.pp.neighbors`) before calling `sc.tl.dpt`.  
- Keep pseudotime metadata standardized: `monocle3_pseudotime`, `sling_pseudotime`, or `dpt_pseudotime`.  
- Compare multiple methods (Monocle3, Slingshot, DPT) for consistency across pseudotime estimates.  
- For large datasets, subsample cells for graph learning to reduce computation time.  
- Record R/Python package versions (`sessionInfo()`, `pip freeze`, `conda env export`) for reproducibility.

---

## 🧾 Attribution
Created and maintained by **Farouk Saaidia (2025)**.  
For research or educational reuse, please cite:

> Saaidia F. (2025). *Single-Cell RNA-seq Trajectory and Pseudotime Analysis Module.*

