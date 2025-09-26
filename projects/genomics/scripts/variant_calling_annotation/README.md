# Variant Calling & Annotation Scripts

## Included Scripts

### Variant Calling
- call_variants_gatk.sh
- call_variants_bcftools.sh
- call_freebayes.sh

### Annotation
- annotate_vep.sh
- annotate_annovar.sh
- annotate_bcftools_csq.sh
- annotate_variants_snpeff.sh

### Post-annotation Utilities
- extract_clinical_variants.sh → Extract variants with clinical significance (ClinVar/COSMIC)
- annotate_population_frequency.sh → Add population allele frequencies from gnomAD or 1000G
- annotate_conservation_scores.sh → Add conservation scores (GERP/PhyloP)

### Pipeline Helpers
- run_all_callers.sh → Wrapper to call all variant callers on the same BAM
- run_all_annotations.sh → Wrapper to sequentially annotate a VCF with all annotation tools
- vcf_cleanup.sh → Cleanup low-quality or duplicate entries from merged VCFs

## Notes
- All scripts take input files as arguments and generate output files.
- Variant caller choice depends on the use-case; annotations add complementary information.
