echo "# FASTQ QC Scripts

| Script Name              | Description                                                | Input                          | Output                        | Notes / Tool Options             |
|--------------------------|------------------------------------------------------------|--------------------------------|-------------------------------|---------------------------------|
| fastqc.sh                | Run FastQC for initial quality assessment                 | R1/R2 FASTQ                    | FastQC HTML and ZIP reports   | Standard FastQC reports         |
| run_multiqc.sh           | Aggregate multiple FastQC reports into a single report    | Folder of FastQC outputs       | MultiQC HTML report           | Provides summary across samples|
| run_trim_galore.sh       | Trim adapters and low-quality bases using Trim Galore      | R1/R2 FASTQ                    | Trimmed FASTQ pair             | Default Trim Galore parameters |
| run_cutadapt.sh          | Adapter and quality trimming using Cutadapt               | R1/R2 FASTQ                    | Trimmed FASTQ pair             | Specify adapters and min length|
| run_trimmomatic.sh       | Adapter trimming, quality filtering, paired/unpaired      | R1/R2 FASTQ                    | Paired & unpaired FASTQ       | Standard Trimmomatic params    |
| run_fastp.sh             | Fast all-in-one trimming and QC                             | R1/R2 FASTQ                    | Trimmed FASTQ pair + HTML report | Fastp default reports          |
| run_all_fastq_qc.sh      | Batch-run any of the above tools on all samples           | Tool name, input dir, output dir | One folder per sample with results | Wrapper for automation        |

**Usage examples:**

- Run FastQC on a single sample:
./fastqc.sh sample_R1.fastq.gz sample_R2.fastq.gz

- Aggregate all FastQC reports with MultiQC:
./run_multiqc.sh ~/data/fastq_qc_results

- Run trimming for all samples with Trim Galore:
./run_all_fastq_qc.sh trim_galore ~/data/fastq ~/data/fastq_qc_results

- Run batch trimming for multiple tools:
./run_all_fastq_qc.sh fastp ~/data/fastq ~/data/fastq_qc_results
" > README.md
