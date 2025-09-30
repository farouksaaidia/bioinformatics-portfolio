# Bulk RNA-seq — Quality Control (QC)

This folder contains scripts to perform **quality control** for bulk RNA-seq (raw FASTQ QC, RNA-seq specific QC on BAMs, library complexity estimation, and aggregation of QC results).  
Each script can operate on a **single file** (FASTQ or BAM) or **all files in a directory** (batch mode).

---

## Scripts Overview

| Script Name | Description | Input | Output | Notes / Tool Options | Example |
|-------------|-------------|-------|--------|----------------------|---------|
| **run_fastqc.sh** | Run FastQC on raw FASTQ files (single or folder) | FASTQ file or directory of FASTQ files | `*_fastqc.html`, `*_fastqc.zip` | Standard per-base/adapter QC; supports gzipped FASTQ | `./run_fastqc.sh /data/fastq /out/fastqc` |
| **run_multiqc.sh** | Aggregate FastQC (and other) reports with MultiQC | Directory containing tool output dirs (FastQC, Picard, RSeQC) | `multiqc_report.html` | Produces a single summary dashboard across samples | `./run_multiqc.sh /out/fastqc /out/multiqc` |
| **run_rseqc.sh** | Run RSeQC RNA-seq checks (e.g., geneBody_coverage, junction saturation) | BAM file or directory of BAMs; reference GTF | RSeQC metric files / plots per sample | Requires sorted & indexed BAMs; RNA-seq specific checks | `./run_rseqc.sh /data/bams ref.gtf /out/rseqc` |
| **run_picard_rnaseq_metrics.sh** | Picard `CollectRnaSeqMetrics` for RNA metrics | BAM, REF_FLAT file, ribosomal intervals | `*_rnaseq_metrics.txt`, summary metrics | Provides rRNA%, coding bases, etc.; needs Picard and REF_FLAT | `./run_picard_rnaseq_metrics.sh sample.bam refFlat.txt ribo.interval /out/picard` |
| **run_preseq.sh** | Estimate library complexity / saturation with Preseq | BAM file or directory of BAMs | `*_preseq.txt` (complexity curves) | Helps decide if more sequencing is useful; preseq installed required | `./run_preseq.sh /data/bams /out/preseq` |
| **run_all_qc.sh** | Wrapper that runs the whole QC pipeline (FastQC→RSeQC→Picard→Preseq→MultiQC) | FASTQ dir, BAM dir, REF_GTF, REF_FLAT, RIBO_INTERVALS, outdir | Combined QC outputs + `multiqc_report.html` | Runs per-sample or batch depending on inputs; aggregates with MultiQC | `./run_all_qc.sh /data/fastq /data/bams ref.gtf refFlat.txt ribo.interval /out/all_qc` |

---

## Usage Notes & Best Practices

- **File naming**: Use a consistent sample naming scheme (e.g. `SAMPLE_R1.fastq.gz` / `SAMPLE_R2.fastq.gz`, or `SAMPLE.bam`) so batch scripts detect pairs automatically.  
- **Paired-end FASTQ**: `run_fastqc.sh` will operate per FASTQ file; pair-aware downstream steps expect `_R1/_R2` naming.  
- **BAM prerequisites**: RSeQC and Picard require **coordinate-sorted and indexed** BAMs (`.bai`). Ensure you run `samtools sort` and `samtools index` or equivalent beforehand.  
- **Reference files**:
  - `RSeQC` uses a GTF for gene body coverage and junction checks.  
  - `CollectRnaSeqMetrics` requires a `refFlat` file and optional ribosomal intervals.  
- **Compute resources**: Increase threads where supported (e.g., FastQC has limited threading; Preseq and Picard can use more resources).  
- **MultiQC**: Always run `run_multiqc.sh` at the end — it provides a single interactive HTML summarizing all QC outputs across samples.  
- **Containers / Reproducibility**: For consistency, run these scripts inside a Conda environment or Docker container with pinned tool versions. Consider adding an `environment.yml` or `Dockerfile` to the project root.  
- **Logs & traceability**: Save stdout/stderr to log files when running batch jobs (e.g., `./run_all_qc.sh ... 2>&1 | tee all_qc.log`) to improve reproducibility and debugging.  
- **Output organization**: Keep QC outputs in a dedicated directory (e.g., `/results/qc/<sample>/`) so downstream scripts can find them easily.

---

## Quick example workflow

```bash
# 1. Raw FASTQ QC
./run_fastqc.sh /data/fastq /results/qc/fastqc

# 2. Align reads and produce sorted/indexed BAMs (see ../bulk_alignment/)
# 3. Run RNA-seq specific QC
./run_rseqc.sh /results/bams ref.gtf /results/qc/rseqc
./run_picard_rnaseq_metrics.sh /results/bams refFlat.txt ribo.interval /results/qc/picard

# 4. Estimate library complexity
./run_preseq.sh /results/bams /results/qc/preseq

# 5. Aggregate with MultiQC
./run_multiqc.sh /results/qc /results/qc/multiqc
