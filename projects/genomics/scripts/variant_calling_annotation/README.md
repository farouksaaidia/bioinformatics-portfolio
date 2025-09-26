# Variant Calling & Annotation Pipeline

This folder contains all scripts related to DNA **variant calling**, **annotation**, and **VCF post-processing**. It is designed for genomic DNA data (BAM files) and provides a complete workflow from raw variant calling to fully annotated VCFs.

---

## 1️⃣ Variant Calling Scripts

| Script | Description | Input | Output | When to Use |
|--------|-------------|-------|--------|------------|
| `call_variants_gatk.sh` | Calls variants using GATK HaplotypeCaller | Sorted BAM, Reference FASTA | Raw VCF | Gold standard for high-quality variant calling; widely used in clinical and research pipelines |
| `call_variants_bcftools.sh` | Calls variants using BCFtools mpileup/call | Sorted BAM, Reference FASTA | Raw VCF | Lightweight alternative; faster for smaller datasets |
| `call_freebayes.sh` | Calls variants using FreeBayes | Sorted BAM, Reference FASTA | Raw VCF | Useful for population-level variant calling and small indels |

**Note:** Each caller generates its own VCF. Use the `run_all_callers.sh` helper to produce individual VCFs and a merged, non-redundant holistic VCF.

---

## 2️⃣ Annotation Scripts

| Script | Description | Input | Output | Notes |
|--------|-------------|-------|--------|------|
| `annotate_vep.sh` | Adds functional annotations using Ensembl VEP | Raw or filtered VCF | Annotated VCF | Adds consequence, gene, and impact information |
| `annotate_annovar.sh` | Adds annotations with ANNOVAR databases | VCF | Annotated VCF | Includes ClinVar, dbSNP, exonic function annotations |
| `annotate_bcftools_csq.sh` | Consequence annotation with BCFtools csq | VCF, Reference FASTA | Annotated VCF | Lightweight annotation of coding consequences |
| `annotate_variants_snpeff.sh` | Predicts variant effects using SnpEff | VCF | Annotated VCF | Functional impact prediction |
| `annotate_population_frequency.sh` | Adds population allele frequencies (gnomAD/1000G) | VCF | Annotated VCF | Helps identify rare vs common variants |
| `annotate_conservation_scores.sh` | Adds conservation metrics (GERP, PhyloP) | VCF | Annotated VCF | Useful to prioritize conserved, potentially functional variants |
| `extract_clinical_variants.sh` | Extracts clinically significant variants | VCF | Filtered VCF | Focus on variants in ClinVar or COSMIC databases |

**Note:** Use `run_all_annotations.sh` to apply all annotators sequentially on a single VCF, producing a fully annotated VCF file.

---

## 3️⃣ Utility / Support Scripts

| Script | Description | Input | Output |
|--------|-------------|-------|--------|
| `filter_vcf.sh` | Filters raw VCFs by quality, depth, or allele frequency | VCF | Filtered VCF |
| `merge_vcfs.sh` | Merges multiple VCFs from different callers or samples | Multiple VCFs | Merged VCF |
| `split_vcf_by_chrom.sh` | Splits large VCFs by chromosome for easier downstream analysis | VCF | Chromosome-specific VCFs |
| `vcf_cleanup.sh` | Cleans up duplicates or low-quality entries after merging | VCF | Cleaned VCF |

---

## 4️⃣ Pipeline Helpers

| Script | Description | Input | Output | User Scenario |
|--------|-------------|-------|--------|--------------|
| `run_all_callers.sh` | Runs all variant callers on a BAM, then merges results | BAM, Reference FASTA | Individual VCFs + merged holistic VCF | Ideal for comprehensive variant detection and cross-validation |
| `run_all_annotations.sh` | Annotates a VCF sequentially using all annotation scripts | Raw or filtered VCF | Fully annotated VCF | Produces a single VCF with all functional, population, and clinical annotations |

---

## Usage Notes

- All scripts expect **sorted BAM files** as input for variant calling.
- Outputs are generally **VCF files**, optionally annotated or filtered.
- Users can pick specific callers or annotators, or use the `run_all_*` helpers for complete workflows.
- Utility scripts help manage large datasets, merge results, and ensure clean, consistent VCFs.
