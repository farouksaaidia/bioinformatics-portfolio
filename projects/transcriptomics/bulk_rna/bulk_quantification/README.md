# Bulk RNA-seq Quantification

## Introduction
Quantification is a critical step in bulk RNA-seq analysis, where the abundance of transcripts or genes is estimated from aligned or raw sequencing reads.  
There are two main strategies:

1. **Alignment-based quantification** – Requires BAM files from alignment tools (e.g., STAR, HISAT2). These methods count reads overlapping annotated features.  
2. **Alignment-free quantification** – Works directly from FASTQ reads, using pseudo-alignment or probabilistic models for faster and more memory-efficient quantification.

---

## Alignment-based Quantification

| Script | Input | Output | Particularity | Use case |
|--------|-------|--------|---------------|----------|
| `run_featurecounts.sh` | BAM file(s), GTF/GFF annotation | Gene-level count matrix | Very fast, widely used; optimized for large datasets | Standard gene count quantification from alignments |
| `run_htseq_count.sh` | BAM file(s), GTF annotation | Gene-level count table | Flexible with different counting modes (union, intersection) | When fine control over counting rules is needed |

---

## Alignment-free Quantification

| Script | Input | Output | Particularity | Use case |
|--------|-------|--------|---------------|----------|
| `run_salmon.sh` | FASTQ file(s), transcriptome index | Transcript-level quantification (TPM, counts) | Fast, bias-aware (GC bias correction) | For transcript-level analysis, fast pipelines |
| `run_kallisto.sh` | FASTQ file(s), transcriptome index | Transcript-level quantification (TPM, counts) | Extremely lightweight, pseudo-alignment | When speed and low memory are critical |
| `run_rsem.sh` | FASTQ file(s), reference (genome or transcriptome) | Transcript & gene-level estimates | Very accurate, benchmarked in many studies | When high accuracy is needed despite heavier compute cost |

---

## Notes
- Use **alignment-based methods** if you already have BAM files and primarily need gene-level counts.  
- Use **alignment-free methods** if you want transcript-level quantification or need faster turnaround.  
- All scripts here have options to run on **single or multiple samples**.  
- Outputs can be directly used for downstream **normalization, differential expression, and visualization**.  
