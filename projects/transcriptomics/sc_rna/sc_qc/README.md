# Single-cell RNA-seq Quality Control (sc_qc)

This folder contains scripts for performing quality control (QC) on single-cell RNA-seq datasets.  
QC steps include initial read quality assessment, empty droplet detection, doublet detection, preprocessing QC checks, and generation of comprehensive QC reports. All scripts support **single-sample or batch-mode workflows** with robust error handling.

| Script Name                  | Input                          | Output                        | When to Use                              | Particularity / Why Use This Script                       |
|-------------------------------|--------------------------------|-------------------------------|-----------------------------------------|----------------------------------------------------------|
| run_fastqc_sc.sh              | R1/R2 FASTQ files              | FastQC HTML & ZIP reports     | Initial read QC                          | Handles single or paired-end reads, standard FastQC     |
| run_multiqc_sc.sh             | Folder of FastQC outputs       | MultiQC HTML report           | Summarize QC across multiple samples    | Provides comprehensive summary across all samples       |
| run_preprocess_seurat.R       | Count matrix / Seurat object   | Filtered / normalized object  | Preprocessing prior to downstream analysis | Uses Seurat workflows, includes filtering and normalization |
| run_preprocess_scanpy.py      | Count matrix / AnnData object  | Filtered / normalized object  | Preprocessing for Scanpy pipelines       | Supports Scanpy workflows, batch mode included          |
| run_emptydrops.R              | Count matrix                   | Filtered matrix               | Remove empty droplets                     | Uses EmptyDrops method for droplet-level filtering       |
| run_scrublet_scanpy.py        | Count matrix                   | Doublet predictions           | Detect doublets in single-cell data      | Single or batch mode, integrates with Scanpy            |
| compare_qc_thresholds.R       | Count matrix + QC metrics      | Heatmaps / plots              | Assess QC thresholds                     | Helps set optimal cutoff values for filtering           |
| generate_qc_report.py         | All QC outputs                 | HTML report                   | Final report generation                  | Aggregates all QC metrics into a professional HTML report |

**Notes:**  
- Always perform QC **before normalization, clustering, or downstream analyses** to remove low-quality cells and technical artifacts.  
- Scripts are designed to work on **single samples or multiple samples** in batch mode.  
- Outputs are standardized to ensure compatibility with Seurat and Scanpy pipelines.  
- This folder focuses on single-cell RNA-seq QC; mitochondrial/ribosomal content or other specialized QC are handled in dedicated folders.

