# sp_visualization

Visualization module for spatial transcriptomics: static publication figures, interactive exploration, report generation, and export of publication-ready assets.

---

## 📂 Folder Purpose

`projects/transcriptomics/spatial_rna/sp_visualization/`  
Produce high-quality static plots, interactive Shiny exploration tools, combined PDF/HTML visual reports, and exportable visual assets (PNGs, embeddings, coordinates) for figure assembly.

---

## 🧩 Scripts Overview

| Script | Language | Input | Output | Description |
|--------|-----------|--------|---------|-------------|
| **plot_spatial_features.R** | R | Seurat `.rds` | PNGs + optional combined PDF | Create spatial feature overlays for genes/metadata (publication-ready) |
| **plot_spatial_domains.R** | R | Seurat `.rds` | PNGs | Plot spatial domains/clusters with summary cluster counts |
| **interactive_spatial_viewer.R** | R (Shiny) | Seurat `.rds` | Interactive app | Lightweight Shiny app for interactive exploration (launch via Rscript) |
| **generate_spatial_visual_report.R** | R | Seurat `.rds` | PDF report | Render an Rmarkdown report summarizing embeddings, domains, and highlighted features |
| **export_visual_assets.sh** | Bash + R | Seurat `.rds` | PNGs + CSVs | Export images, embeddings, and coordinate CSVs for figure panels |
| **.gitkeep** | text | — | — | Keeps folder tracked by Git |

---

## ⚙️ Example Workflows

### 1️⃣ Make spatial feature plots (PNG + combined PDF)
Rscript projects/transcriptomics/spatial_rna/sp_visualization/plot_spatial_features.R \
  -i results/normalized/SAMPLE01_sct.rds \
  -f "GeneA,GeneB,SPP1" \
  -o results/visuals/SAMPLE01_features \
  --png --pdf

### 2️⃣ Plot domains and cluster counts
Rscript projects/transcriptomics/spatial_rna/sp_visualization/plot_spatial_domains.R \
  -i results/clustering/SAMPLE01_clusters.rds \
  -c bayes_domains \
  -o results/visuals/SAMPLE01_domains

### 3️⃣ Launch interactive viewer
Rscript projects/transcriptomics/spatial_rna/sp_visualization/interactive_spatial_viewer.R \
  results/preprocessing/SAMPLE01_seurat.rds

### 4️⃣ Generate combined visual report
Rscript projects/transcriptomics/spatial_rna/sp_visualization/generate_spatial_visual_report.R \
  -i results/preprocessing/SAMPLE01_seurat.rds \
  -o reports/SAMPLE01_visual_report.pdf \
  -f "GENE1,GENE2"

### 5️⃣ Export visual assets for figure assembly
projects/transcriptomics/spatial_rna/sp_visualization/export_visual_assets.sh \
  -i results/preprocessing/SAMPLE01_seurat.rds -o results/visual_assets/SAMPLE01

---

## 🧠 Best Practices

- Use **high-resolution PNGs** (300 dpi) for publication figures; keep separate vector exports (PDF/SVG) where possible.  
- When overlaying features on histology, visually inspect alignment and cropping.  
- Use the interactive viewer to confirm gene/cluster spatial patterns before exporting static figures.  
- For multi-panel figures, export consistent color scales and a shared legend.  
- Save intermediate assets (embeddings, coords) to reproduce figures outside R (e.g., Illustrator, Python).  
- Document figure-generation steps in a small script so figures are reproducible.

---

## 🧾 Attribution

Created and maintained by **Farouk Saaidia (2025)**.
