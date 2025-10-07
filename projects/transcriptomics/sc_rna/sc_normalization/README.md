# Single-Cell RNA-seq Normalization

This folder contains scripts to normalize single-cell RNA-seq data for both single and multiple samples. It includes classical log-normalization, SCTransform, Scran-based normalization, scaling/regression, multi-sample normalization wrappers, and batch correction methods.

Normalization is a critical step to remove technical variation and prepare data for downstream analysis such as clustering, dimensionality reduction, and differential expression.

---

## Single-Sample Normalization Scripts

| Script | Input | Output | When to Use | Particularity |
|--------|-------|--------|-------------|---------------|
| run_log_normalization.R | Seurat RDS object | Normalized Seurat RDS object | Standard log-normalization | Simple, widely used, fast |
| run_sctransform.R | Seurat RDS object | SCTransformed Seurat RDS | High variance genes or batch variation | Models technical noise explicitly |
| run_scran_pooling.R | SingleCellExperiment RDS | Normalized SCE object | Complex, sparse datasets | Pooling-based, handles library size differences |
| run_scaling.R | Seurat RDS object | Scaled Seurat RDS object | Before PCA/clustering | Regress out unwanted covariates |

---

## Multi-Sample / Batch Normalization Scripts

| Script | Input | Output | When to Use | Particularity |
|--------|-------|--------|-------------|---------------|
| run_multi_sample_normalization.R | Folder of RDS files | Folder of normalized RDS files | Multiple samples, same method | Wrapper for log or SCTransform |
| run_batch_correction.R | Folder of RDS files | Combined RDS object | Multiple batches | Supports Harmony or Seurat integration |

---

### Notes
- Always inspect QC metrics before and after normalization.
- Multi-sample normalization can handle different methods (`log` or `SCTransform`) to harmonize processing across datasets.
- Batch correction should be applied **after normalization** but **before clustering/dimensionality reduction**.
- Error handling is included in all scripts to ensure reproducibility and robust pipeline integration.

