# Short-Read RNA-seq: Alignment and Quantification

## Introduction
This module performs genome alignment or transcript quantification for short-read RNA-seq data using STAR, HISAT2, or Salmon. It produces BAM files or quant folders required for downstream normalization and differential expression.

## Scripts

| Script | Description | Input | Output | Command | When to Use |
|--------|-------------|--------|---------|----------|--------------|
| run_star_alignment.sh | Aligns paired-end reads using STAR and outputs sorted BAM files. | FASTQ files, STAR index, sample sheet | Sorted + indexed BAM files | bash run_star_alignment.sh -r ref -s samples.tsv -o out -t 8 | Standard genome alignment with splice junction detection. |
| run_hisat2_alignment.sh | Aligns paired-end reads using HISAT2 and sorts BAM with samtools. | FASTQ files, HISAT2 index, sample sheet | Sorted + indexed BAM files | bash run_hisat2_alignment.sh -x index -s samples.tsv -o out -t 8 | Alternative splice-aware aligner; lighter memory usage. |
| run_salmon_quant.sh | Performs alignment-free transcript quantification. | FASTQ files, Salmon index, sample sheet | Salmon quant directories | bash run_salmon_quant.sh -i idx -s samples.tsv -o out -t 8 | When transcript-level quantification or speed is required. |
| merge_bams.sh | Merges multiple BAMs into one file. | Directory of BAMs | Merged BAM + index | bash merge_bams.sh -i bam_dir -o merged.bam | For combined QC or visualization across samples. |

