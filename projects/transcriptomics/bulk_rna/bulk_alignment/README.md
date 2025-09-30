# Bulk RNA-seq Alignment Scripts

This folder contains scripts for aligning bulk RNA-seq FASTQ samples to a reference genome.
It includes single-sample alignment scripts for STAR, HISAT2, and Bowtie2, as well as a batch wrapper to align multiple samples in one go.
These scripts produce sorted and indexed BAM files suitable for downstream quantification and analysis.

| Script Name                | Description                                | Input                                      | Output                        | Notes / Tool Options                                       |
|----------------------------|--------------------------------------------|-------------------------------------------|-------------------------------|------------------------------------------------------------|
| align_star.sh              | Align a single sample using STAR           | FASTQ R1/R2, STAR genome index, out dir  | Sorted & indexed BAM          | High-speed, splice-aware aligner, ideal for most bulk RNA-seq datasets |
| align_hisat2.sh            | Align a single sample using HISAT2         | FASTQ R1/R2, HISAT2 genome index, out dir| Sorted & indexed BAM          | Splice-aware aligner, good for large genomes and exon junction mapping |
| align_bowtie2.sh           | Align a single sample using Bowtie2        | FASTQ R1/R2, Bowtie2 genome index, out dir| Sorted & indexed BAM         | Less sensitive to splicing, faster for unspliced or targeted RNA |
| run_all_bulk_alignment.sh  | Batch-run all alignment tools on multiple samples | List of FASTQ pairs, reference genome indexes, output folder | Sorted & indexed BAM per sample/tool | Automates multi-sample alignment, allows tool selection |

**Usage Examples:**

- Align a single sample with STAR:
  ./align_star.sh sample_R1.fastq.gz sample_R2.fastq.gz /path/to/star_index /output/dir

- Align a single sample with HISAT2:
  ./align_hisat2.sh sample_R1.fastq.gz sample_R2.fastq.gz /path/to/hisat2_index /output/dir

- Align multiple samples using the batch wrapper:
  ./run_all_bulk_alignment.sh sample_list.txt /path/to/indexes /output/dir STAR,HISAT2

**Notes:**

- Ensure the reference genome indexes exist before running any alignment script.
- Adjust the number of threads (-p) as needed based on your system resources.
- Output BAM files are sorted and indexed automatically for downstream tools like featureCounts or StringTie.

