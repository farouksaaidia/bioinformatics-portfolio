#!/bin/bash
# run_all_bulk_alignment.sh → Batch-align multiple bulk RNA-seq samples
# Input: alignment tool (star/hisat2/bowtie2), input folder, reference, output folder
# Output: Sorted BAMs for all samples

TOOL=$1
INPUT_DIR=$2
REF=$3
OUT_DIR=$4

mkdir -p $OUT_DIR

for SAMPLE in $INPUT_DIR/*_R1*.fastq*; do
    BASE=$(basename $SAMPLE _R1.fastq.gz)
    FASTQ1="$INPUT_DIR/${BASE}_R1.fastq.gz"
    FASTQ2="$INPUT_DIR/${BASE}_R2.fastq.gz"
    
    case $TOOL in
        star)
            bash bulk_rna/bulk_alignment/align_star.sh $FASTQ1 $FASTQ2 $REF $OUT_DIR/$BASE
            ;;
        hisat2)
            bash bulk_rna/bulk_alignment/align_hisat2.sh $FASTQ1 $FASTQ2 $REF $OUT_DIR/$BASE
            ;;
        bowtie2)
            bash bulk_rna/bulk_alignment/align_bowtie2.sh $FASTQ1 $FASTQ2 $REF $OUT_DIR/$BASE
            ;;
        *)
            echo "Unsupported tool: $TOOL"
            exit 1
            ;;
    esac
done
