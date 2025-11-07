# 🧬 Short-Read RNA-seq Preprocessing

## Overview
This module performs the **initial stage of short-read RNA-seq analysis**, preparing high-quality, adapter-free FASTQ files for downstream alignment and quantification.  
It ensures sequencing integrity, removes adapters, and validates sample metadata before mapping.

---

## 🔬 Pipeline Stage: Preprocessing

**Purpose:**  
Ensure all input FASTQ files are high-quality and properly formatted for alignment.

**Key Tasks:**
1. Perform quality control using **FastQC**.  
2. Aggregate QC results with **MultiQC**.  
3. Trim adapters and low-quality reads using **Cutadapt**.  
4. Validate the sample sheet to confirm data integrity.

---

## ⚙️ Scripts and Utilities

| Script | Description | Input | Output | When to Use |
|--------|--------------|--------|---------|--------------|
| run_fastqc_sr.sh | Runs FastQC on all FASTQ files to generate quality metrics. | Raw FASTQ files | FastQC HTML and ZIP reports | Immediately after receiving raw FASTQs. |
| run_multiqc_sr.sh | Aggregates FastQC results into one report. | FastQC reports directory | multiqc_report.html | After running FastQC to view a summary. |
| run_trimming_sr.sh | Trims adapters and low-quality bases (SE/PE). | FASTQ files + adapter file | Trimmed FASTQs | After QC, before alignment. |
| validate_sample_sheet.py | Validates that all FASTQs in the sample sheet exist. | Sample sheet TSV | Validation log | Before alignment to avoid missing samples. |

---

## 🧠 Example Commands

Run FastQC  
bash run_fastqc_sr.sh -i data/fastq_raw -o results/qc/fastqc -t 8

Run MultiQC  
bash run_multiqc_sr.sh -i results/qc/fastqc -o results/qc/multiqc

Run Trimming (Single-End)  
bash run_trimming_sr.sh -i data/fastq_raw -o results/trimmed -a adapters.fa -m SE -t 8

Run Trimming (Paired-End)  
bash run_trimming_sr.sh -i data/fastq_raw -o results/trimmed -a adapters.fa -m PE -t 8

Validate Sample Sheet  
python3 validate_sample_sheet.py -i metadata/sample_sheet.tsv

---

## 📁 Outputs

| File or Folder | Description |
|----------------|--------------|
| fastqc/ | Individual FastQC reports |
| multiqc_report.html | Combined QC summary |
| trimmed/ | Cleaned FASTQs ready for mapping |
| validated_samples.log | Log of validated samples |

---

## 🧩 Dependencies

| Tool | Purpose |
|------|----------|
| FastQC | Read quality assessment |
| MultiQC | Combine QC outputs |
| Cutadapt | Adapter and quality trimming |
| Python 3 + pandas | Sample validation |

---

