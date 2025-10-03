# Bulk RNA-seq Analysis Pipeline

This folder contains scripts and workflows for **Bulk RNA-seq analysis**, structured step-by-step to reflect a professional transcriptomics pipeline.  
Each subfolder corresponds to one stage of analysis, containing alternative approaches and tools commonly used in bioinformatics.

---

## 📌 Pipeline Overview

1. **Alignment** → Map RNA-seq reads to a reference genome/transcriptome  
2. **Quality Control (QC)** → Assess raw and aligned data quality  
3. **Post-alignment QC** → Inspect alignment metrics, biases, and coverage  
4. **Quantification** → Count reads per gene or transcript  
5. **Normalization** → Adjust for sequencing depth, library size, composition bias  
6. **Differential Expression** → Identify differentially expressed genes/transcripts  
7. **Visualization** → Generate publication-ready plots  

---

## 📂 Subfolder Summaries

### 1. Alignment (`bulk_alignment/`)

| Script | Input | Output | When to Use | Particularity |
|--------|--------|---------|-------------|---------------|
| `align_hisat2.sh` | FASTQ, reference index | SAM/BAM | Splice-aware alignment | Popular for RNA-seq, handles splice junctions |
| `align_star.sh` | FASTQ, reference index | SAM/BAM | High-performance RNA-seq alignment | Extremely fast, sensitive |
| `align_salmon.sh` | FASTQ, transcriptome index | Quant files | Alignment-free | Direct quantification, lightweight |
| `align_kallisto.sh` | FASTQ, transcriptome index | Quant files | Alignment-free | Pseudoalignment, very fast |
| `align_bowtie2.sh` | FASTQ, reference index | SAM/BAM | General-purpose | Older, less RNA-specific |
| `align_bwa.sh` | FASTQ, reference index | SAM/BAM | Rare in RNA-seq | Mostly for DNA, shown for completeness |

---

### 2. Quality Control (`bulk_qc/`)

| Script | Input | Output | When to Use | Particularity |
|--------|--------|---------|-------------|---------------|
| `run_fastqc.sh` | FASTQ | HTML report | Single-sample QC | Standard per-sample QC |
| `run_multiqc.sh` | QC reports | Combined HTML | Summarize QC across samples | Aggregates FastQC, etc. |
| `run_fastp.sh` | FASTQ | Cleaned FASTQ, QC HTML | QC + trimming | All-in-one solution |
| `run_trimmomatic.sh` | FASTQ | Trimmed FASTQ | Flexible trimming | Widely used |
| `run_cutadapt.sh` | FASTQ | Trimmed FASTQ | Adapter removal | Highly customizable |
| `run_trim_galore.sh` | FASTQ | Trimmed FASTQ | Adapter + QC | Wrapper for Cutadapt + FastQC |

---

### 3. Post-alignment QC (`bulk_post_alignment_qc/`)

| Script | Input | Output | When to Use | Particularity |
|--------|--------|---------|-------------|---------------|
| `run_samtools_flagstat.sh` | BAM | Stats file | Quick QC metrics | Lightweight |
| `run_samtools_stats.sh` | BAM | Detailed stats | Deeper QC | Rich stats |
| `run_bamqc.sh` | BAM | HTML report | Visual QC | Part of Qualimap |
| `run_rseqc_geneBody.sh` | BAM, GTF | Plots | Coverage bias | Detects 5’/3’ bias |
| `run_rseqc_inner_distance.sh` | BAM, GTF | Plots | Insert size distribution | Library QC |
| `run_rseqc_infer_experiment.sh` | BAM, GTF | Orientation info | Strand-specific QC | Ensures correct strandedness |
| `run_picard_markdups.sh` | BAM | Metrics + deduped BAM | Duplicate detection | Library complexity check |

---

### 4. Quantification (`bulk_quantification/`)

| Script | Input | Output | When to Use | Particularity |
|--------|--------|---------|-------------|---------------|
| `run_featurecounts.sh` | BAM, GTF | Count matrix | Gene-level | Fast + accurate |
| `run_htseq.sh` | BAM, GTF | Count matrix | Gene-level | Flexible |
| `run_salmon.sh` | FASTQ, index | Quant files | Transcript-level | Alignment-free |
| `run_kallisto.sh` | FASTQ, index | Quant files | Transcript-level | Very fast |

---

### 5. Normalization (`bulk_normalization/`)

| Script | Input | Output | When to Use | Particularity |
|--------|--------|---------|-------------|---------------|
| `run_tpm.sh` | Count matrix | Normalized matrix | Transcript-level normalization | Per-million scaling |
| `run_fpkm.sh` | Count matrix | Normalized matrix | Gene/transcript expression | Accounts for gene length |
| `run_cpm.sh` | Count matrix | Normalized matrix | Basic normalization | Library-size correction |
| `run_deseq2_norm.R` | Count matrix | Normalized matrix | DESeq2 workflows | Variance-stabilized |
| `run_edger_tmm.R` | Count matrix | Normalized matrix | edgeR workflows | Trimmed mean normalization |
| `run_voom_norm.R` | Count matrix | Normalized matrix | limma workflows | Converts counts → log2 |
| `run_upperquartile_norm.R` | Count matrix | Normalized matrix | Robust normalization | Reduces skew |
| `run_quantile_norm.R` | Count matrix | Normalized matrix | Distribution normalization | Matches distributions |
| `run_rlog_norm.R` | Count matrix | Normalized matrix | Small-sample datasets | Stabilizes variance |

---

### 6. Differential Expression (`bulk_differential_expression/`)

| Script | Input | Output | When to Use | Particularity |
|--------|--------|---------|-------------|---------------|
| `run_deseq2.R` | Counts, metadata | DEG table | Standard | Most widely used |
| `run_edger.R` | Counts, metadata | DEG table | Flexible | Handles replicates |
| `run_limma_voom.R` | Counts, metadata | DEG table | Microarray-like | Voom transform |
| `run_ballgown.R` | Transcript quant | DEG table | Transcript-level | Isoform analysis |
| `run_sleuth.R` | Kallisto quant | DEG table | Alignment-free | Lightweight |

**Advanced methods:**

| Script | Input | Output | When to Use | Particularity |
|--------|--------|---------|-------------|---------------|
| `run_dream.R` | Counts, covariates | DEG table | Mixed models | Complex designs |
| `run_mait.R` | Counts | DEG table | Multi-omics | Integrates data |
| `run_wgcna.R` | Expression matrix | Modules | Network biology | Gene co-expression |
| `run_rankprod.R` | Counts | DEG table | Small sample sizes | Non-parametric |
| `run_bayseq.R` | Counts | DEG table | Bayesian analysis | Alternative modeling |

---

### 7. Visualization (`bulk_visualization/`)

| Script | Input | Output | When to Use | Particularity |
|--------|--------|---------|-------------|---------------|
| `run_pca_plot.R` | Normalized matrix | PCA plot | Sample clustering | Detect batch effects |
| `run_heatmap.R` | Normalized matrix | Heatmap | Expression overview | Top DEGs |
| `run_volcano_plot.R` | DEG table | Volcano plot | DE visualization | Effect size vs significance |
| `run_ma_plot.R` | DEG table | MA plot | Expression shift | DE overview |
| `run_sample_distance_heatmap.R` | Normalized matrix | Heatmap | Sample relationships | Outlier detection |
| `run_barplot_topgenes.R` | DEG table | Barplot | Top DEGs | Publication-ready |
| `run_pathway_enrichment_plot.R` | Enrichment results | Plot | Pathway visualization | Functional analysis |
| `run_gsea_plot.R` | GSEA results | Plot | Gene set enrichment | Pathway-level |
| `run_boxplot_genes.R` | Counts, metadata | Boxplots | Expression per condition | Validate DE |
| `run_upset_plot.R` | Multiple DEG tables | UpSet plot | DEG overlap | Multi-condition comparison |

---

## 📌 General Notes

- Each folder includes **multiple alternative tools** to showcase expertise.  
- Scripts are structured for **single sample or multiple sample batch processing**.  
- The pipeline is modular → recruiters can see how I can adapt workflows depending on needs.  
- Where multiple tools exist, their **particular strengths and tradeoffs** are highlighted.  

