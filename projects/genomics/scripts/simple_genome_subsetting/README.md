# Simple Genome Subsetting

This module contains a collection of scripts for extracting specific subsets of genomic data.  
They are useful when working with large sequencing files and you only need to focus on particular regions, samples, or variants without processing the entire dataset.  

---

## Scripts Overview

| Script Name              | Function (What it does)                                  | Input                          | Output                        | Why / When to Use                                                                 | Example Usage |
|--------------------------|----------------------------------------------------------|--------------------------------|-------------------------------|-----------------------------------------------------------------------------------|---------------|
| extract_chromosome.sh    | Extracts a specific chromosome from a FASTA/FASTQ/BAM/VCF | FASTA/FASTQ/BAM/VCF            | Subset file with only one chromosome | When analyzing chromosome-specific signals or reducing file size.                   | `./extract_chromosome.sh input.bam chr1 chr1_subset.bam` |
| extract_region.sh        | Extracts a genomic region by coordinates                 | BAM/VCF + region (chr:start-end)| Subset BAM/VCF file          | To focus on candidate regions (e.g., GWAS hits, exons).                            | `./extract_region.sh input.vcf chr1:1000-2000 region_subset.vcf` |
| extract_random_reads.sh  | Randomly selects a subset of reads                       | FASTQ                          | Subsampled FASTQ              | Useful for testing workflows on smaller data, QC checks, or downsampling.          | `./extract_random_reads.sh input.fastq 10000 subset.fastq` |
| subset_reads_by_name.sh  | Extracts reads matching a list of read names             | FASTQ/BAM + list of read names | Subset FASTQ/BAM              | To isolate reads of interest (e.g., specific barcodes, outliers).                  | `./subset_reads_by_name.sh input.bam names.txt subset.bam` |
| subset_variants_vcf.sh   | Subsets variants from a VCF based on criteria (region, ID, type) | VCF + filtering parameters   | Filtered VCF                 | When you only want specific variants (e.g., SNPs, INDELs, candidate variants).     | `./subset_variants_vcf.sh input.vcf snps_only.vcf --snps-only` |
| subset_samples_vcf.sh    | Extracts specific samples from a multi-sample VCF        | Multi-sample VCF + sample list | Subset VCF with selected samples | To focus on individuals of interest in population or cohort studies.               | `./subset_samples_vcf.sh cohort.vcf sample_list.txt subset.vcf` |

---

## Notes

- All scripts assume **bgzip/tabix indexing** where applicable (`.vcf.gz`, `.bam`).
- Outputs are compatible with downstream tools (samtools, bcftools, bedtools, etc.).
- Designed to be **lightweight helpers** for exploratory analysis and pipeline preparation.

