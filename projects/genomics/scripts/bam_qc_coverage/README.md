# BAM QC and Coverage Analysis

This folder contains scripts for **post-alignment quality control (QC)** and **coverage analysis**.  
These steps help evaluate mapping quality, sequencing depth, and uniformity before downstream variant calling.  

---

## 📊 Scripts Overview

| Script Name                  | Description | Input | Output | When to Use / Why |
|------------------------------|-------------|-------|--------|-------------------|
| run_picard_collectmetrics.sh | Collect alignment summary metrics using **Picard** (insert size, % mapped reads, duplication rate). | Sorted BAM, Reference FASTA | Metrics TXT, optional PDF | Use when you need detailed alignment QC beyond basic stats. |
| run_samtools_stats.sh        | Generate alignment statistics with **samtools stats**. | Sorted BAM | Stats TXT | Lightweight QC overview; run quickly on all BAMs. |
| run_qc_multiqc.sh            | Aggregate multiple QC reports (Picard + Samtools) into a **single MultiQC report**. | Folder of QC outputs | MultiQC HTML report | Use to summarize results across many samples. |
| run_samtools_depth.sh        | Compute per-base sequencing depth across genome/exome. | BAM, Reference FASTA | Depth TXT | Use for whole-genome coverage visualization and depth distribution. |
| run_bedtools_coverage.sh     | Calculate coverage over specific regions (e.g., exome capture kit targets). | BAM, BED file | Coverage TXT | Use when working with targeted sequencing or panels. |
| run_mosdepth.sh              | Fast genome-wide or region-based coverage calculation with **mosdepth**. | BAM, Reference FASTA | Depth summaries (TXT, BED, bedgraph) | Faster and more memory-efficient than samtools depth; use for large WGS. |
| run_cov_qc_summary.sh        | Combine coverage outputs into summary tables (% bases >10x, mean depth, etc.). | Depth results from other tools | Summary TXT/CSV | Use to produce consolidated coverage QC metrics for downstream reporting. |
| run_all_bam_qc_coverage.sh   | Wrapper to run **all BAM QC + coverage scripts** in one go. | BAM files, Reference FASTA, BED (optional) | Aggregated QC + coverage reports | Use for batch automation across large cohorts. |
