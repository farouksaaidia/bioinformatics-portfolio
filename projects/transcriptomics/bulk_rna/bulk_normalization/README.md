# Bulk RNA-seq Normalization

## Introduction
This folder contains scripts for performing various normalization methods on bulk RNA-seq data.  
Normalization is critical to account for sequencing depth, gene length, and technical variability.  
We include both **scaling-based methods** (CPM, TPM, FPKM) and **statistical/model-based methods** (DESeq2, edgeR, limma, quantile normalization, VST, rlog).  

---

## Scaling-based Normalization

| Script            | Input                  | Output                 | Particularity                          | Use Case |
|-------------------|-----------------------|-----------------------|----------------------------------------|----------|
| `run_cpm.sh`      | Raw counts matrix     | CPM normalized counts | Divides counts by library size, scaled | Quick normalization, simple comparisons |
| `run_tpm.sh`      | Counts + gene lengths | TPM values            | Accounts for gene length and depth     | Within-sample expression comparison |
| `run_fpkm.sh`     | Counts + gene lengths | FPKM values           | Similar to TPM, older method           | Historically used, mostly replaced |

---

## Model/statistical Normalization

| Script                  | Input                           | Output                  | Particularity                                 | Use Case |
|--------------------------|--------------------------------|------------------------|-----------------------------------------------|----------|
| `run_deseq2_norm.R`      | Counts + sample info           | Normalized counts      | Median ratio normalization (robust)           | Downstream DE analysis |
| `run_edger_tmm.R`        | Counts + sample info           | Normalized counts      | TMM normalization (robust scaling)            | DE with edgeR |
| `run_voom_norm.R`        | Counts + sample info           | logCPM + weights       | Prepares data for limma linear models         | Linear modeling, DE |
| `run_quantile_norm.R`    | Counts matrix                  | Quantile-normalized    | Forces identical distributions across samples | Cross-dataset comparisons |
| `run_vst.R`              | Counts + sample info           | Variance stabilized    | Reduces variance dependency on mean           | PCA, clustering, viz |
| `run_rlog.R`             | Counts + sample info           | rlog-transformed data  | Regularized log, best for small datasets      | Visualization, PCA |

---

## Notes
- **Scaling methods (CPM/TPM/FPKM)** are simple, fast, but not suitable for statistical testing.  
- **Model/statistical methods (DESeq2, TMM, voom, VST, rlog)** are recommended for **differential expression** and robust downstream analyses.  
- **Quantile normalization** is powerful for harmonizing datasets, but may remove biological signal if overused.  

