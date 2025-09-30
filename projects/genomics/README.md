# Genomics Pipeline Overview

This repository contains a comprehensive genomics analysis pipeline, structured for professional bioinformatics workflows. The folder structure is modular, separating tasks for clarity and reproducibility.

---

## Folder and Script Overview

| Subfolder                   | Description |
|------------------------------|-------------|
| `data`                       | Contains raw FASTQ files, reference genomes, and analysis results. |
| `scripts/alignment`          | Scripts for sequence alignment using BWA, Bowtie2, Minimap2, and BLAST. |
| `scripts/fastq_qc`           | Initial QC and trimming of FASTQ files using FastQC, MultiQC, Trim Galore, Cutadapt, Trimmomatic, and fastp. |
| `scripts/bam_qc_coverage`    | BAM QC and coverage analysis scripts using Samtools, Picard, Qualimap, and Bedtools. |
| `scripts/file_conversion_compression` | Scripts for file format conversions, compression, and indexing (SAM/BAM/VCF/FASTQ/FASTA). |
| `scripts/reference_indexing` | Scripts to index reference genomes for various alignment tools (BWA, Bowtie2, Picard, Samtools). |
| `scripts/simple_genome_subsetting` | Scripts to subset genome sequences, reads, or variants based on regions, chromosomes, or samples. |
| `scripts/variant_calling_annotation` | Variant calling (GATK, FreeBayes, BCFtools) and annotation pipelines (VEP, SnpEff, Annovar, BCFtools). Includes utility scripts for merging, filtering, and post-processing. |
| `scripts/visualization_prep` | Scripts to prepare BAM/VCF files for visualization, convert to BigWig/BedGraph, and subset regions. |

---

## General Notes

- All scripts are designed to be **modular and reusable**.
- Input and output directories are clearly defined; outputs are compatible with downstream analyses.
- Utility scripts provide automation and batch processing capabilities.
- Each subfolder contains its own `README.md` with detailed explanations of scripts, usage, inputs, outputs, and user scenarios.

---

