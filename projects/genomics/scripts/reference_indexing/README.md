# Reference Indexing

This folder contains scripts to generate reference genome indexes required for DNA sequence alignment and downstream analyses.  
Indexing is a critical preprocessing step that allows aligners and tools to efficiently map sequencing reads or access genomic regions.

---

## 📑 Scripts Overview

| Script Name          | Tool / Function      | Input                  | Output                              | When to Use                                                                 | Usage Example |
|----------------------|----------------------|------------------------|-------------------------------------|----------------------------------------------------------------------------|---------------|
| bwa_index.sh         | BWA indexing         | Reference FASTA        | `.amb`, `.ann`, `.bwt`, `.pac`, `.sa` | Use when preparing a reference for **BWA alignment** of DNA reads.           | `./bwa_index.sh reference.fasta` |
| bowtie2_index.sh     | Bowtie2 indexing     | Reference FASTA        | `.bt2` index files                   | Use when preparing a reference for **Bowtie2 alignment** of DNA reads.       | `./bowtie2_index.sh reference.fasta index_prefix` |
| samtools_faidx.sh    | Samtools FASTA index | Reference FASTA        | `.fai` index file                    | Use when needing **random access** to reference sequences (samtools, bcftools). | `./samtools_faidx.sh reference.fasta` |
| picard_create_dict.sh| Picard dictionary    | Reference FASTA        | `.dict` file                         | Required for **GATK/Picard pipelines** (variant calling, processing).        | `./picard_create_dict.sh reference.fasta` |
| index_cleanup.sh     | Cleanup script       | Old index files        | Cleaned directory                    | Run before regenerating indexes to avoid conflicts.                         | `./index_cleanup.sh reference.fasta` |

---

## 📝 Notes

- Keep **FASTA, `.fai`, `.dict`, and aligner-specific indexes** in the same directory for consistency.  
- `bwa` and `bowtie2` indexes are **not interchangeable**.  
- Use `index_cleanup.sh` to avoid mixing old and new index files.  
- RNA-seq aligners (e.g., STAR, HISAT2) belong to the **transcriptomics** folder, not here.  

