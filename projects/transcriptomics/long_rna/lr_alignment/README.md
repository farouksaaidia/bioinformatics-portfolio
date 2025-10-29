# Long-Read RNA-seq — Alignment

This module performs splice-aware mapping of long RNA reads to a reference genome.
It supports multiple samples, automated indexing, and alignment QC reporting.

Goals:
- Align full-length RNA reads to the genome with long-read optimized methods
- Sort and index alignments for downstream isoform analysis
- Compute mapping statistics and genome coverage metrics
- Support high-throughput automation via sample manifests

-------------------------------------
SPOTLIGHT ON LONG-READ ALIGNMENT

Long reads span full transcripts and introns.
Therefore, mapping needs:
- Splice-aware algorithms
- Soft-clipping tolerance (long read errors)
- Chimeric detection for structural RNA events

Minimap2 is highly recommended for this purpose.

-------------------------------------
AVAILABLE SCRIPTS

| Script | Function | Input | Output | Use When |
|--------|----------|-------|--------|----------|
| run_minimap2_align.sh | Splice-aware mapping (Minimap2) | FASTQ.gz or sample list | Sorted BAM + BAI index | Standard single- or multi-sample alignment |
| sort_index_bam.sh | Sorting + indexing only | BAM | Sorted BAM + BAI | Cleanup or after external alignment |
| multi_sample_manifest.sh | Batch alignment automation | TSV manifest | Multiple BAM files | Many samples to align consistently |
| generate_alignment_qc.R | Mapping statistics report | BAM | CSV + PDF report | Validate alignment performance |
| extract_alignment_logs.sh | Parse mapping logs | Alignment log | Text summary | Quick metrics for dashboards |

-------------------------------------
EXAMPLE WORKFLOWS

Single sample:
./run_minimap2_align.sh -r reference.fa -i reads.fastq.gz -o aligned/

Batch mode:
manifest.tsv format:
reads1.fastq.gz   reference.fa
reads2.fastq.gz   reference.fa

./multi_sample_manifest.sh manifest.tsv aligned/

QC reporting:
./generate_alignment_qc.R aligned/sample.bam aligned/qc/

-------------------------------------
OUTPUT ARTIFACTS

Aligned BAM files (sorted + indexed)
Mapping statistics (CSV + PDF)
Log summaries with mapped percentage

-------------------------------------
NEXT PIPELINE STEP

Proceed to:
lr_transcript_reconstruction
(isoform detection and transcriptome annotation)

-------------------------------------
Maintainer: Farouk Saaidia
Last updated: 2025-10
