# Short-Read RNA-seq: Normalization

## Introduction
This module normalizes raw or quantified gene counts from short-read RNA-seq. It supports DESeq2 size-factor normalization, edgeR CPM normalization, TPM conversion, and merging Salmon quantification outputs.

## Scripts

| Script | Description | Input | Output | Command | When to Use |
|--------|-------------|--------|---------|----------|--------------|
| normalize_counts_deseq2.R | Normalizes counts using DESeq2 size factors. | Raw counts CSV/TSV + metadata TSV | normalized_counts_deseq2.csv | Rscript normalize_counts_deseq2.R -c counts.csv -m meta.tsv -o out | When performing DESeq2-based differential expression. |
| normalize_counts_edger_cpm.R | Generates CPM-normalized counts using edgeR. | Raw counts CSV/TSV | normalized_counts_cpm.csv | Rscript normalize_counts_edger_cpm.R -c counts.csv -o out | When quick CPM normalization is needed for clustering or QC. |
| counts_to_tpm_sr.R | Converts raw counts to TPM using gene lengths. | Count matrix + gene length file | counts_tpm.csv | Rscript counts_to_tpm_sr.R -c counts.csv -l lengths.tsv -o out | When comparing gene expression across samples. |
| merge_salmon_counts.py | Merges Salmon quant.sf files into matrices. | Salmon quant directories | salmon_counts_matrix.csv, salmon_tpm_matrix.csv | python3 merge_salmon_counts.py -i quant_dir -o out | When using Salmon quantification instead of alignment-based counts. |

