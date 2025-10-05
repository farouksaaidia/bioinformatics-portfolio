# 🧬 Single-cell RNA-seq Preprocessing

This folder contains scripts for the **preprocessing** of single-cell RNA-seq data — transforming raw FASTQ files into clean, usable gene-count matrices.

---

## 📋 Table: Preprocessing Scripts Overview

| Script | Input | Output | When to Use | Particularity |
|--------|--------|---------|--------------|----------------|
| **run_cellranger_mkfastq.sh** | BCL files (raw sequencer output) | FASTQ files | First step if you have Illumina BCLs | Converts raw BCLs to FASTQ using Cell Ranger |
| **run_cellranger_count.sh** | FASTQ files + reference genome | Filtered gene-cell matrix | Standard 10x Genomics pipeline | Combines alignment, filtering, and UMI counting |
| **run_fastqc_sc.sh** | FASTQ files | HTML reports | Always first QC check | Checks read quality per cell/sample |
| **run_multiqc_sc.sh** | FastQC reports | Combined HTML report | After QC | Aggregates QC from multiple samples |
| **run_emptydrops.R** | Gene-cell count matrix | Filtered matrix with valid cells | After CellRanger / Alevin | Detects empty droplets in droplet-based data |
| **run_preprocess_scanpy.py** | Count matrix (H5AD/MTX) | Filtered + normalized Scanpy object | For Python users | Scanpy-based preprocessing (filtering, normalization, log1p) |
| **run_preprocess_seurat.R** | Count matrix (CSV/MTX) | Filtered Seurat object | For R users | Seurat-based preprocessing (filtering, normalization, variable gene detection) |

---

## 🧠 Notes
- Always perform **FASTQ-level QC** (FastQC + MultiQC) before preprocessing.  
- Choose **CellRanger**, **STARsolo**, or **Alevin** depending on your dataset source and protocol.  
- **EmptyDrops** helps exclude background barcodes from droplet-based single-cell experiments.  
- Preprocessing results feed directly into **normalization** and **integration** workflows downstream.

