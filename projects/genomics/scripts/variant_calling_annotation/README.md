# Variant Calling & Annotation Scripts

This folder contains scripts for DNA variant calling, annotation, and VCF post-processing.

## Variant Calling
- call_variants_gatk.sh → GATK HaplotypeCaller
- call_variants_bcftools.sh → bcftools mpileup + call
- call_freebayes.sh → FreeBayes
- run_all_callers.sh → Run all callers and produce separate + merged non-redundant VCF

## Annotation
- annotate_vep.sh → Ensembl VEP
- annotate_annovar.sh → ANNOVAR
- annotate_variants_snpeff.sh → SnpEff
- annotate_bcftools_csq.sh → bcftools CSQ
- annotate_population_frequency.sh → Add gnomAD / 1000G allele frequencies
- extract_clinical_variants.sh → ClinVar / COSMIC extraction
- annotate_conservation_scores.sh → Add GERP / PhyloP scores
- run_all_annotations.sh → Run all annotators sequentially to produce one fully annotated VCF

## Post-processing
- vcf_cleanup.sh → Remove duplicates / low-quality entries

## Usage Examples
\`\`\`bash
# Run all callers
./run_all_callers.sh sample.bam reference.fa output_dir

# Annotate a VCF
./run_all_annotations.sh output_dir/merged.vcf.gz output_dir/annotated

# Clean a VCF
./vcf_cleanup.sh output_dir/annotated/annotated_final.vcf output_dir/annotated/annotated_final_cleaned.vcf
\`\`\`
