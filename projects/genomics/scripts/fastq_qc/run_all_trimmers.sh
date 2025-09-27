#!/bin/bash
# Universal FASTQ QC batch runner
# USAGE: ./run_all_fastq_qc.sh <tool> <input_dir> <output_dir>
# TOOL: trim_galore / cutadapt / trimmomatic / fastp

TOOL=$1
INPUT_DIR=$2
OUT_DIR=$3

for fq1 in $INPUT_DIR/*_R1.fastq.gz; do
    fq2=${fq1/_R1/_R2}
    sample=$(basename $fq1 _R1.fastq.gz)
    mkdir -p $OUT_DIR/$sample

    case $TOOL in
        trim_galore)
            ~/bioinformatics-portfolio/projects/genomics/scripts/fastq_qc/run_trim_galore.sh $fq1 $fq2 $OUT_DIR/$sample
            ;;
        cutadapt)
            OUT1=$OUT_DIR/$sample/${sample}_R1_trimmed.fastq.gz
            OUT2=$OUT_DIR/$sample/${sample}_R2_trimmed.fastq.gz
            ~/bioinformatics-portfolio/projects/genomics/scripts/fastq_qc/run_cutadapt.sh $fq1 $fq2 $OUT1 $OUT2
            ;;
        trimmomatic)
            OUT1P=$OUT_DIR/$sample/${sample}_R1_paired.fastq.gz
            OUT1U=$OUT_DIR/$sample/${sample}_R1_unpaired.fastq.gz
            OUT2P=$OUT_DIR/$sample/${sample}_R2_paired.fastq.gz
            OUT2U=$OUT_DIR/$sample/${sample}_R2_unpaired.fastq.gz
            ~/bioinformatics-portfolio/projects/genomics/scripts/fastq_qc/run_trimmomatic.sh $fq1 $fq2 $OUT1P $OUT1U $OUT2P $OUT2U
            ;;
        fastp)
            OUT1=$OUT_DIR/$sample/${sample}_R1_trimmed.fastq.gz
            OUT2=$OUT_DIR/$sample/${sample}_R2_trimmed.fastq.gz
            REPORT=$OUT_DIR/$sample/${sample}_fastp_report.html
            ~/bioinformatics-portfolio/projects/genomics/scripts/fastq_qc/run_fastp.sh $fq1 $fq2 $OUT1 $OUT2 $REPORT
            ;;
        *)
            echo "Unknown tool: $TOOL"
            exit 1
            ;;
    esac
done
