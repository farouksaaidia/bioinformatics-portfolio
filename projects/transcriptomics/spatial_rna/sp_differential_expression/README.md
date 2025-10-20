# sp_differential_expression

Differential expression and spatial enrichment analysis for spatial transcriptomics datasets (Visium, Slide-seq, Stereo-seq).  
This module detects marker genes per spatial domain, performs spatial autocorrelation testing (Moran's I), computes pairwise differential expression between spatial domains, runs pathway enrichment for DEG sets, and compiles comprehensive reports.

---

## 📂 Folder Purpose

`projects/transcriptomics/spatial_rna/sp_differential_expression/`  
Contains scripts to compute domain/cluster marker genes, test spatial autocorrelation of gene expression, run pairwise domain DE analyses, perform pathway enrichment, and generate summarized PDF reports.

---

## 🧩 Scripts Overview

| Script | Language | Input | Output | Description |
|--------|----------|-------|--------|-------------|
| **find_spatial_markers.R** | R | clustered Seurat `.rds` | per-cluster marker CSVs | Identify marker genes per spatial cluster/domain using Seurat `FindMarkers` |
| **spatial_moran_test.R** | R | Seurat `.rds` | CSV of Moran's I stats | Test spatial autocorrelation for genes using Moran's I |
| **differential_between_domains.R** | R | Seurat `.rds` | pairwise DE CSVs | Pairwise differential expression comparisons between spatial domains |
| **pathway_enrichment_spatial.R** | R | DEG CSV directory | GO/KEGG CSVs | Run functional enrichment on DEG lists using `clusterProfiler` |
| **generate_spatial_deg_report.R** | R | Seurat `.rds` + DEG dir | PDF report | Compile heatmaps, volcanoes, and enrichment summaries into a PDF |
| **.gitkeep** | text | — | — | Keeps folder tracked by Git |

---

## ⚙️ Example Workflows

### 1️⃣ Find spatial markers per domain
Rscript projects/transcriptomics/spatial_rna/sp_differential_expression/find_spatial_markers.R \
  -i results/clustering/SAMPLE01_clusters.rds \
  -c bayes_domains \
  -o results/differential/SAMPLE01_markers

### 2️⃣ Moran's I spatial autocorrelation
Rscript projects/transcriptomics/spatial_rna/sp_differential_expression/spatial_moran_test.R \
  -i results/normalized/SAMPLE01_sct.rds \
  -g "GENE1,GENE2" \
  -k 6 \
  -o results/differential/SAMPLE01_moran.csv

### 3️⃣ Pairwise differential expression between domains
Rscript projects/transcriptomics/spatial_rna/sp_differential_expression/differential_between_domains.R \
  -i results/clustering/SAMPLE01_domains.rds \
  -d bayes_domains \
  -o results/differential/SAMPLE01_pairwise_DE

### 4️⃣ Pathway enrichment for DEGs
Rscript projects/transcriptomics/spatial_rna/sp_differential_expression/pathway_enrichment_spatial.R \
  -d results/differential/SAMPLE01_pairwise_DE \
  -o results/differential/SAMPLE01_enrichment

### 5️⃣ Generate combined DEG & enrichment report
Rscript projects/transcriptomics/spatial_rna/sp_differential_expression/generate_spatial_deg_report.R \
  -s results/clustering/SAMPLE01_clusters.rds \
  -d results/differential/SAMPLE01_pairwise_DE \
  -e results/differential/SAMPLE01_enrichment \
  -o reports/SAMPLE01_spatial_DEG_report.pdf

---

## 🧠 Best Practices

- Use normalized and scaled assay (SCT or scran) when running DE tests.  
- Apply multiple-testing correction (`p_val_adj`) and sensible thresholds (e.g., `FDR < 0.05`, `|log2FC| > 0.25`).  
- Moran's I is computationally expensive; select a subset of genes or the top variable genes to test.  
- Validate domain DE results against histology and spatial overlays.  
- For enrichment, provide a background gene list (expressed genes) to reduce false positives.  
- Record all package versions (`sessionInfo()`) for reproducibility.

---

## 🧾 Attribution

Created and maintained by **Farouk Saaidia (2025)**.
