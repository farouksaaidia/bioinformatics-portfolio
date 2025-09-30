# Visualization Preparation Scripts

This folder contains scripts to prepare genomic data files (**BAM, VCF**) for visualization in tools like **IGV**, **UCSC Genome Browser**, and **JBrowse**.  
These scripts ensure files are **sorted, indexed, and converted** into formats that visualization tools can efficiently read.

---

## 🧬 BAM Visualization Prep

| Script Name             | Input            | Output            | Purpose & Use Case | Visualization Tool |
|--------------------------|------------------|-------------------|--------------------|--------------------|
| **sort_bam.sh**         | BAM              | Sorted BAM        | Sorts BAM by coordinates to enable indexing and efficient browsing. | IGV, JBrowse |
| **index_bam.sh**        | Sorted BAM       | BAM index (.bai)  | Creates BAM index required to load BAMs in genome browsers. | IGV, UCSC, JBrowse |
| **subset_bam_region.sh**| BAM + region     | Subset BAM        | Extracts reads from a specific genomic region (e.g., gene locus). | IGV (zoomed regions), UCSC |
| **bam_to_bedgraph.sh**  | BAM              | bedGraph          | Converts BAM to bedGraph format to show coverage profiles. | UCSC Genome Browser |
| **bedgraph_to_bigwig.sh**| bedGraph + chrom.sizes | BigWig    | Converts bedGraph to BigWig (compressed, efficient coverage format). | UCSC, IGV, JBrowse |

---

## 🧬 VCF Visualization Prep

| Script Name             | Input     | Output              | Purpose & Use Case | Visualization Tool |
|--------------------------|-----------|---------------------|--------------------|--------------------|
| **index_vcf.sh**        | VCF       | bgzipped VCF + .tbi | Compresses and indexes VCF files so they can be queried and visualized. | IGV, JBrowse, UCSC |

---

## 📌 Usage Notes

- **BAM Sorting & Indexing**  
  - Always sort BAM files before indexing.  
  - Required for **IGV** and most genome browsers.  

- **Coverage Tracks**  
  - Convert BAM → bedGraph → BigWig to generate **coverage plots**.  
  - BigWig is preferred because it’s smaller and faster than bedGraph.  
  - Used in **UCSC Genome Browser**, **IGV**, and **JBrowse**.  

- **VCF Indexing**  
  - Compress and index VCFs with `bgzip` + `tabix`.  
  - Allows browsers to **fetch variants by region** without loading the full file.  
  - Supported by **IGV**, **JBrowse**, and **UCSC Genome Browser**.  

