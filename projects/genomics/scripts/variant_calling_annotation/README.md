# Variant Calling & Annotation Scripts

This folder contains scripts for **DNA variant calling and annotation**.

## Included Scripts

### Variant Calling
- **call_gatk.sh** - GATK HaplotypeCaller
- **call_freebayes.sh** - FreeBayes
- **call_bcftools.sh** - Samtools mpileup + BCFtools

### Annotation
- **annotator1.sh**, **annotator2.sh**, etc. (each script adds specific annotation fields to VCF)

### Support / Utility
- **filter_vcf.sh** → Filter raw VCFs by quality, depth, allele frequency.
- **merge_vcfs.sh** → Merge multiple VCF files from different callers/samples.
- **split_vcf_by_chrom.sh** → Split large VCFs per chromosome for downstream analysis.

## Inputs & Outputs
- **Variant Calling Scripts:** Input = sorted BAM, reference FASTA; Output = raw VCF.
- **Annotation Scripts:** Input = VCF; Output = annotated VCF.
- **Support Scripts:** Input = VCF(s); Output = filtered, merged, or split VCF(s).

## Usage Recommendations
- **Variant caller choice:** Depends on project goals (sensitivity, specific variant types).
- **Annotation:** Each annotator adds complementary data (allele frequency, gene impact, clinical relevance).
- **Support scripts:** Used to clean, merge, or split VCFs before/after annotation.

