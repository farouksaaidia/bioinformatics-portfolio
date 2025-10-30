# 📊 Long RNA-seq Normalization

## Overview
This module handles normalization and batch correction for bulk or long-read RNA-seq datasets. It converts raw count matrices into comparable expression units (TPM, CPM, RPKM, or TMM) and corrects technical variation across samples or batches.

Normalization ensures that observed differences in gene expression reflect biological variation rather than sequencing depth, transcript length, or batch effects.

---

## 🔬 Pipeline Stage: Normalization
Purpose: Transform raw counts into comparable scales and adjust for systematic biases.  
This step typically follows alignment and precedes differential expression analysis.

Key tasks include:
1. Converting raw counts to TPM/CPM/RPKM for within-sample normalization.  
2. Computing TMM factors for cross-sample scaling.  
3. Applying batch correction (ComBat-seq).  
4. Generating diagnostic visualizations of normalization quality.

---

## ⚙️ Scripts and Utilities

| Script | Description | Input | Output | When to Use |
|--------|--------------|--------|---------|--------------|
| counts_to_tpm.R | Converts transcript-level counts to TPM, CPM, or RPKM. Handles optional length files. | Raw transcript counts (TSV) + optional lengths | Normalized matrix (TPM/CPM/RPKM) | Use for transcript-level data (especially long-read sequencing like Iso-Seq). |
| compute_normalizations.R | Computes gene-level TPM, CPM, RPKM, and TMM normalization. | Counts + optional gene lengths | Normalized matrices (counts_tpm.csv, counts_cpm_tmm.csv, etc.) | Use for standard short/long-read RNA-seq normalization. |
| combat_seq_batch_correction.R | Applies ComBat-seq (or limma fallback) to correct batch effects in count data. | Raw counts + metadata (sample, batch) | Batch-corrected matrix (counts_combatseq_corrected.csv) | Use when samples come from different sequencing runs or experimental batches. |
| normalize_batch_wrapper.sh | Wrapper to normalize multiple samples using a manifest or a single matrix. | Manifest or counts matrix + optional gene lengths | Organized output per sample or combined | Use for high-throughput or automated normalization pipelines. |
| visualize_normalization.R | Generates PCA, boxplots, and density plots before and after normalization. | Raw and normalized matrices | PDF reports (boxplots_pre_post_log2.pdf, pca_pre_post.pdf, etc.) | Use for QC after normalization or batch correction. |

---

## 🧠 Example Commands

Example 1: TPM normalization (transcript-level)
Rscript counts_to_tpm.R --counts data/transcript_counts.tsv --lengths data/transcript_lengths.tsv --output results/normalization --method TPM

Example 2: Gene-level normalization and TMM scaling
Rscript compute_normalizations.R --counts data/gene_counts.tsv --lengths data/gene_lengths.tsv --methods tpm,cpm,tmm --outdir results/normalization

Example 3: Batch correction (ComBat-seq)
Rscript combat_seq_batch_correction.R --counts results/normalization/counts_tpm.csv --metadata data/sample_metadata.tsv --outdir results/batch_corrected

Example 4: Multi-sample wrapper
bash normalize_batch_wrapper.sh -c data/gene_counts.tsv -l data/gene_lengths.tsv -o results/normalization -m "tpm,cpm,tmm"

Example 5: Visualization of normalization results
Rscript visualize_normalization.R --raw data/gene_counts.tsv --normalized results/normalization/counts_tpm.csv --outdir results/plots

---

## 📁 Outputs
All normalized and corrected matrices are saved under the chosen --outdir, including:
- counts_tpm.csv  
- counts_cpm.csv  
- counts_rpkm.csv  
- tmm_factors.csv  
- counts_combatseq_corrected.csv  
- Diagnostic plots (PDFs)

---

## 🧩 Dependencies
- R packages: optparse, dplyr, readr, edgeR, sva, limma, ggplot2, ggrepel  
- Runtime: R ≥ 4.0.0  
- Optional: manifest TSV for batch mode normalization

---

## 🧭 Next Step
Proceed to the next stage: lr_differential_expression, where normalized expression data are used to identify significantly differentially expressed genes.
