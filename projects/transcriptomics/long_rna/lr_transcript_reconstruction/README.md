# Long-Read RNA-seq — Transcript Reconstruction and Quantification

This module reconstructs and quantifies full-length transcripts from long-read RNA-seq data.
It integrates multiple specialized tools to:
- Assemble full isoforms (not fragmented exons)
- Identify novel splicing events and gene models
- Quantify transcript-level expression
- Assess structural and quality metrics

-------------------------------------
OVERVIEW

Long-read sequencing (ONT, PacBio) captures entire transcripts.
After alignment, these reads can be used to **build isoform-level annotations** and measure expression accurately.

Typical workflow:
FLAIR → StringTie2 → TALON → SQANTI3 → Unified Isoform Report

1. **FLAIR**: Correct and collapse reads into transcripts.
2. **StringTie2 (long mode)**: Quantify isoforms per gene.
3. **TALON**: Build and manage an annotated transcript database.
4. **SQANTI3**: Perform QC, classify novel isoforms, and remove artifacts.
5. **Isoform report**: Generate summary statistics and visualizations.

-------------------------------------
AVAILABLE SCRIPTS

| Script | Function | Input | Output | When to Use |
|--------|----------|-------|--------|--------------|
| run_flair_pipeline.sh | Isoform correction and reconstruction using FLAIR | BAM (aligned reads), genome.fa, annotation.gtf | Collapsed isoforms (GTF) | When starting from raw aligned reads; ideal for identifying new transcripts and splicing events |
| run_stringtie2_long.sh | Quantify transcript expression using StringTie2 (long-read mode) | BAM, annotation.gtf | Gene abundance table (TSV) and GTF | When you want expression values per isoform after reconstruction |
| run_talon_annotation.sh | Build a TALON transcript database and quantify known + novel isoforms | BAM, reference GTF, TALON DB | Annotated database + abundance files | When you want consistent annotation tracking across experiments |
| run_sqanti3_qc.sh | Quality control and structural classification using SQANTI3 | Assembled GTF, genome.fa, reference GTF | QC metrics, classification tables, filtered GTF | When you want to validate and filter isoforms before final reporting |
| generate_isoform_report.R | Summarize isoform categories and QC statistics into figures and tables | SQANTI3 classification file | PDF + CSV summary | When you need a final visual summary for publication or review |

-------------------------------------
EXAMPLE WORKFLOW

Example (FLAIR → StringTie2 → SQANTI3):

1. FLAIR:
   ./run_flair_pipeline.sh -g genome.fa -a annotation.gtf -b sample.bam -o flair_out/
2. Quantify with StringTie2:
   ./run_stringtie2_long.sh -b flair_out/isoforms.bam -a flair_out/isoforms.gtf -o stringtie_out/
3. Run SQANTI3 QC:
   ./run_sqanti3_qc.sh flair_out/isoforms.gtf genome.fa annotation.gtf sqanti_out/
4. Summarize results:
   ./generate_isoform_report.R sqanti_out/classification.txt isoform_report/

-------------------------------------
OUTPUT ARTIFACTS

| Output | Description |
|--------|--------------|
| isoforms.gtf | Full-length transcript models |
| gene_abundance.tsv | Isoform-level expression |
| talon.db | Transcript database (known + novel) |
| sqanti3_classification.txt | Isoform quality and structure |
| isoform_types.pdf | Visualization of isoform composition |
| category_counts.csv | Table of isoform structural categories |

-------------------------------------
NOTES

- FLAIR and TALON can both reconstruct transcripts; choose based on your preferred format or reproducibility needs.
- SQANTI3 can process output from either FLAIR or TALON.
- All tools are compatible with ONT and PacBio long reads.
- Outputs are standardized for integration with differential isoform expression and functional enrichment downstream.

-------------------------------------
NEXT PIPELINE STEP

Proceed to:
lr_isoform_expression
(isoform-level expression analysis and visualization)

-------------------------------------
Maintainer: Farouk Saaidia
Last updated: 2025-10
