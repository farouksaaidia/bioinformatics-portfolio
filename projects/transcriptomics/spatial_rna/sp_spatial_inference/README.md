# sp_spatial_inference

Spatial inference module for ligand–receptor communication, neighborhood analysis, and spatial interaction network discovery.  
This step decodes intercellular communication patterns and tissue-level signaling structures from spatial transcriptomics data.

---

## 📂 Folder Purpose

`projects/transcriptomics/spatial_rna/sp_spatial_inference/`  
Implements ligand–receptor inference (NicheNet, LIANA), neighborhood enrichment (Giotto, Squidpy), method comparison, and network visualization.

---

## 🧩 Scripts Overview

| Script | Language | Input | Output | Description |
|--------|-----------|--------|---------|-------------|
| **run_nichenet_inference.R** | R | annotated Seurat `.rds` | ligand activity CSV | Predict ligand–receptor interactions between sender & receiver cell types |
| **run_liana_analysis.py** | Python | `.h5ad` | CSV | Run LIANA consensus ligand–receptor inference |
| **run_giotto_neighborhood.R** | R | Giotto `.rds` | CSV | Compute neighborhood enrichment across cell types |
| **run_squidpy_interactions.py** | Python | `.h5ad` | CSV | Analyze spatial co-occurrence and proximity patterns |
| **compare_inference_methods.py** | Python | multiple CSVs | merged comparison + overlaps | Compare ligand–receptor results from multiple tools |
| **generate_inference_report.R** | R | merged CSV | PDF | Visualize ligand–receptor network and signal flow |
| **.gitkeep** | text | — | — | Keeps folder tracked by Git |

---

## ⚙️ Example Workflows

### 1️⃣ NicheNet inference
Rscript projects/transcriptomics/spatial_rna/sp_spatial_inference/run_nichenet_inference.R \
  -i results/mapping/SAMPLE01_annotated.rds \
  -s "Fibroblast,Endothelial" \
  -r "T_cell" \
  -o results/inference/SAMPLE01_nichenet.csv

### 2️⃣ LIANA consensus inference
projects/transcriptomics/spatial_rna/sp_spatial_inference/run_liana_analysis.py \
  -i results/mapping/SAMPLE01.h5ad \
  -o results/inference/SAMPLE01_liana.csv

### 3️⃣ Giotto neighborhood enrichment
Rscript projects/transcriptomics/spatial_rna/sp_spatial_inference/run_giotto_neighborhood.R \
  -i results/spatial/SAMPLE01_giotto.rds \
  -o results/inference/SAMPLE01_neighborhood.csv

### 4️⃣ Squidpy spatial interactions
projects/transcriptomics/spatial_rna/sp_spatial_inference/run_squidpy_interactions.py \
  -i results/spatial/SAMPLE01.h5ad \
  -o results/inference/SAMPLE01_squidpy

### 5️⃣ Multi-method comparison
projects/transcriptomics/spatial_rna/sp_spatial_inference/compare_inference_methods.py \
  --nichenet results/inference/SAMPLE01_nichenet.csv \
  --liana results/inference/SAMPLE01_liana.csv \
  --giotto results/inference/SAMPLE01_neighborhood.csv \
  --out_prefix results/inference/SAMPLE01_compare

### 6️⃣ Generate network report
Rscript projects/transcriptomics/spatial_rna/sp_spatial_inference/generate_inference_report.R \
  -i results/inference/SAMPLE01_compare_merged.csv \
  -o reports/SAMPLE01_inference_network.pdf

---

## 🧠 Best Practices

- Always ensure cell-type annotations are accurate before running inference.  
- Run **multiple inference tools** (e.g. NicheNet, LIANA, Giotto) — compare overlapping ligand–receptor pairs.  
- Use **pathway-level validation** (from enrichment modules) to confirm biological relevance.  
- Visualize **networks and spatial context together** to identify signaling hubs.  
- For tumor or immune contexts, inspect immune checkpoint and cytokine-related interactions.  
- Keep versioned ligand–receptor databases for reproducibility.

---

## 🧾 Attribution

Created and maintained by **Farouk Saaidia (2025)**.
