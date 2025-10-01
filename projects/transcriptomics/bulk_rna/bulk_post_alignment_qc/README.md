# Bulk RNA-seq Post-alignment QC

This folder contains scripts for performing post-alignment quality control (QC) on bulk RNA-seq BAM files. The QC metrics help ensure data quality before downstream analyses such as quantification, normalization, and differential expression.

| Script Name                     | Description                                                 | Input                  | Output                                    | Notes / When to Use                                         | Usage Example                                               |
|---------------------------------|-------------------------------------------------------------|-----------------------|------------------------------------------|------------------------------------------------------------|------------------------------------------------------------|
| run_bamqc.sh                     | Wrapper to run multiple BAM QC tools sequentially          | BAM file               | Folder with all QC outputs               | Automates running all QC scripts for a BAM                 | ./run_bamqc.sh sample.bam                                    |
| run_picard_markdups.sh           | Mark duplicate reads using Picard                           | BAM file               | BAM file with duplicates marked           | Use after alignment to mark PCR duplicates                | ./run_picard_markdups.sh sample.bam sample_marked.bam      |
| run_rseqc_geneBody.sh             | Gene body coverage analysis using RSeQC                     | BAM file, GTF file     | Gene body coverage plots                  | Evaluate uniformity of coverage across gene body         | ./run_rseqc_geneBody.sh sample.bam annotation.gtf          |
| run_rseqc_infer_experiment.sh     | Infer strandedness of RNA-seq library                       | BAM file, GTF file     | Strandedness report                        | Check if library is stranded or unstranded               | ./run_rseqc_infer_experiment.sh sample.bam annotation.gtf  |
| run_rseqc_inner_distance.sh       | Inner distance distribution between paired reads            | BAM file, GTF file     | Insert size distribution plots            | Useful for checking library prep and fragment sizes       | ./run_rseqc_inner_distance.sh sample.bam annotation.gtf    |
| run_samtools_flagstat.sh          | Quick summary of BAM file alignment statistics             | BAM file               | Text summary report                        | Fast overview of mapped reads, proper pairs, duplicates   | ./run_samtools_flagstat.sh sample.bam                       |
| run_samtools_stats.sh             | Detailed BAM QC metrics                                     | BAM file               | Text stats file + plots with plot-bamstats | Deep QC: coverage, GC content, error rates, insert size  | ./run_samtools_stats.sh sample.bam                           |

**Notes:**
- Scripts are intended to be run after alignment and optional duplicate marking.
- `run_bamqc.sh` can be used as a one-stop wrapper to generate all QC outputs.
- Use `samtools flagstat` for quick checks and `samtools stats` for detailed assessment.
- RSeQC scripts require a proper GTF annotation file.
