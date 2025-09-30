# 📂 File Conversion & Compression Utilities

This folder contains **utility scripts for converting between genomic file formats, compressing files, and indexing them**.  
These tools are essential in every bioinformatics workflow to ensure **data efficiency, interoperability, and accessibility** for downstream analyses.

---

## 📊 Scripts Overview

| Script Name                      | Description                                           | Input                        | Output                         | When to Use / Why                                                                 | Example Usage |
|----------------------------------|-------------------------------------------------------|------------------------------|--------------------------------|-----------------------------------------------------------------------------------|---------------|
| **fastq_to_fasta.sh**            | Convert FASTQ reads to FASTA format                   | FASTQ (R1/R2, gz or plain)   | FASTA                          | When you only need sequences without quality scores (e.g., BLAST, motif search)   | `./fastq_to_fasta.sh sample.fastq.gz sample.fasta` |
| **sam_to_bam.sh**                | Convert SAM to BAM                                    | SAM                          | BAM                            | To reduce size and speed up processing; BAM is compressed and binary               | `./sam_to_bam.sh input.sam output.bam` |
| **bam_to_sam.sh**                | Convert BAM back to SAM                               | BAM                          | SAM                            | When human-readable format is required                                             | `./bam_to_sam.sh input.bam output.sam` |
| **bam_to_fastq.sh**              | Convert BAM to FASTQ                                  | BAM                          | FASTQ R1, FASTQ R2             | To recover FASTQ reads from aligned BAM (e.g., re-alignment, QC)                  | `./bam_to_fastq.sh aln.bam R1.fastq R2.fastq` |
| **compress_bgzip.sh**            | Compress a file with bgzip (block compression)        | VCF, BED, etc.               | .gz (bgzipped file)            | Needed for large genomic files that require random access                         | `./compress_bgzip.sh variants.vcf variants.vcf.gz` |
| **index_tabix.sh**               | Create tabix index for bgzipped files                 | .vcf.gz, .bed.gz, etc.       | .tbi index                     | Enables fast region queries (`bcftools view -r chr1:1000-2000`)                   | `./index_tabix.sh variants.vcf.gz` |
| **convert_vcf_bcf.sh**           | Convert between VCF and BCF formats                   | VCF or BCF                   | BCF or VCF                     | BCF is faster/lighter for large datasets; VCF is human-readable                   | `./convert_vcf_bcf.sh input.vcf output.bcf` |
| **compress_decompress_gzip.sh**  | Compress or decompress with gzip                      | Any file / .gz file          | .gz or decompressed file       | General-purpose compression or extraction                                         | `./compress_decompress_gzip.sh compress file.txt` |
| **bam_sort_index.sh**            | Sort and index a BAM file                             | BAM                          | Sorted BAM + .bai index        | Essential before variant calling, IGV visualization, or downstream processing     | `./bam_sort_index.sh input.bam sorted.bam` |

---

## 🚀 Usage

1. Make scripts executable once (if not already):
   ```bash
   chmod +x *.sh
Run the script with its required arguments (examples in table above).

📝 Notes & Best Practices
Always sort and index BAM files before running downstream analysis (variant calling, IGV).

Use bgzip + tabix for all large genomic files (VCF, BED, GFF) if you plan region-based queries.

Keep both human-readable (SAM/VCF) and compressed/binary (BAM/BCF) versions if your workflow requires debugging + performance.

fastq_to_fasta should only be used when quality scores are irrelevant (otherwise you lose crucial QC info).

Use bam_to_fastq sparingly, as re-generating FASTQ from BAM may not perfectly reproduce the original FASTQ.

