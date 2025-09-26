# Variant Calling & Annotation Scripts

This folder contains scripts for **DNA variant calling and annotation** using real and mock datasets.

## Variant Calling Scripts
1. **call_variants_gatk.sh**  
   - Uses GATK HaplotypeCaller  
   - Input: sorted BAM, reference FASTA  
   - Output: raw VCF  
   - Best for high-quality WGS/WES samples  

2. **call_freebayes.sh**  
   - Uses FreeBayes  
   - Input: sorted BAM, reference FASTA  
   - Output: raw VCF  
   - Works well for pooled or population samples  

3. **call_variants_bcftools.sh**  
   - Samtools mpileup + BCFtools call  
   - Input: sorted BAM, reference FASTA  
   - Output: raw VCF  
   - Fast and lightweight for small datasets  

## Annotation Scripts
1. **annotate_vep.sh**  
   - Annotates VCF using Ensembl VEP  
   - Adds gene impact, consequence, and protein effect  

2. **annotate_annovar.sh**  
   - Annotates VCF with ANNOVAR  
   - Adds population frequencies, clinical significance, and gene function  

3. **annotate_bcftools_csq.sh**  
   - Uses bcftools csq for basic consequences  
   - Adds canonical transcript consequence annotations  

4. **annotate_variants_snpeff.sh**  
   - Annotates variants with SnpEff  
   - Adds predicted effects, impact, and functional class  

## Support / Utility Scripts
1. **filter_vcf.sh**  
   - Filters raw VCFs based on QUAL, DP, AF  

2. **merge_vcfs.sh**  
   - Merges multiple VCFs from different callers/samples  

3. **split_vcf_by_chrom.sh**  
   - Splits large VCFs by chromosome for easier downstream analysis  

## Usage Recommendations
- **Variant caller choice:** Depends on project goals (sensitivity, specific variant types, pooled vs. single sample)  
- **Annotation scripts:** Each adds complementary biological and clinical information  
- **Support scripts:** Used for cleaning, merging, or splitting VCFs before/after annotation  

