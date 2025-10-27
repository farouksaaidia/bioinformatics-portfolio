# Long-Read RNA-seq — Preprocessing

This module prepares raw long-read sequencing data from Oxford Nanopore (ONT) and PacBio for downstream alignment, isoform detection, and quantification.

Goals:
- Convert raw signals into reads
- Remove adapters and artifacts
- Assess read quality and yield
- Extract biological features like Poly(A) tails
- Standardize output into FASTQ or CCS BAM files

-------------------------------------
PLATFORM DIFFERENCES

| Feature | ONT FAST5 | ONT FASTQ | PacBio Subreads | PacBio HiFi FASTQ |
|--------|-----------|-----------|----------------|------------------|
| Basecalling required | Yes | No | No | No |
| Adapter trimming recommended | Yes | Yes | Optional | Optional |
| Poly(A) tail detection | Yes | Limited | Yes | Yes |
| CCS generation required | No | No | Yes | No |

-------------------------------------
AVAILABLE SCRIPTS (THIS FOLDER)

| Script | Platform | Input | Output | Purpose | Use When |
|--------|----------|-------|--------|---------|----------|
| run_guppy_basecalling.sh | ONT | FAST5 folder | FASTQ.gz + logs | GPU-enabled basecalling | You have FAST5 |
| run_porechop.sh | ONT | FASTQ.gz | trimmed FASTQ.gz | Remove adapters + chimera | Before QC or alignment |
| run_nanoplot_qc.sh | ONT or PacBio | FASTQ.gz | QC HTML report + plots | Measure read length/quality | After trimming or CCS |
| run_polya_tail_detection.sh | ONT (FAST5) | FAST5 folder | CSV | Poly(A) tail length estimation | Transcript stability studies |
| run_pacbio_ccs.sh | PacBio | subreads.bam | ccs.bam | Generate high-accuracy consensus reads | Always with PacBio subreads |
| run_ccs_to_fastq.sh | PacBio | ccs.bam | FASTQ.gz | Standardize CCS reads | Before mapping downstream |

-------------------------------------
COMMON WORKFLOWS

ONT workflow:
1. Basecalling → run_guppy_basecalling.sh
2. Trimming → run_porechop.sh
3. QC → run_nanoplot_qc.sh
4. Optional PolyA → run_polya_tail_detection.sh

PacBio workflow:
1. CCS generation → run_pacbio_ccs.sh
2. Convert CCS → run_ccs_to_fastq.sh
3. QC → run_nanoplot_qc.sh

-------------------------------------
OUTPUT SUMMARY

FASTQ.gz → standardized reads for mapping  
CCS BAM → high-accuracy PacBio reads  
QC report (HTML) → read distributions, N50, throughput  
PolyA CSV → estimated Poly(A) tail lengths  

-------------------------------------
NEXT STEP

Proceed to:
lr_alignment
(splice-aware long-read mapping)

-------------------------------------
Maintainer: Farouk Saaidia
Last updated: 2025-10
