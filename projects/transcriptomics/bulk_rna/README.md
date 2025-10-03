# Bulk RNA-seq Differential Expression

This folder contains scripts to perform **differential gene expression (DGE) analysis** from bulk RNA-seq data.  
Different statistical frameworks are implemented to allow flexibility, cross-validation, and demonstration of expertise with widely used methods.

---

## 📌 Core Methods (Most Commonly Used in Practice)

| Script                | Input                                | Output                       | When to Use                                                   | Particularity Compared to Others |
|------------------------|--------------------------------------|------------------------------|---------------------------------------------------------------|----------------------------------|
| **run_deseq2_DE.R**    | Count matrix + sample metadata       | DEGs table (log2FC, p-value) | Standard choice for bulk RNA-seq with replicates               | Robust normalization (size factors), shrinkage for variance |
| **run_deseq2_lfcShrink.R** | Same as above                   | DEGs with shrunken log2FC    | When accurate log fold-changes are needed (e.g. small samples) | Uses shrinkage estimators (apeglm/ashr) for better effect sizes |
| **run_deseq2_multi.R** | Count matrix + metadata (multi-cond) | DEGs across multiple groups  | For complex designs with >2 conditions                        | Handles multi-factor and interactions |
| **run_edger_DE.R**     | Count matrix + sample metadata       | DEGs table                   | Very popular alternative to DESeq2                             | Empirical Bayes dispersion estimation, robust for small replicates |
| **run_limma_voom.R**   | Count matrix + sample metadata       | DEGs table                   | When data has many samples or complex designs                  | Uses precision weights with linear models, very fast and memory-efficient |
| **run_metaDE.R**       | Results from multiple methods        | Consensus DEG table          | When comparing DE methods or integrating multiple outputs       | Provides meta-analysis across DESeq2, edgeR, limma, etc. |

---

## 📌 Advanced / Specialized Methods

| Script                | Input                                | Output                       | When to Use                                                   | Particularity Compared to Others |
|------------------------|--------------------------------------|------------------------------|---------------------------------------------------------------|----------------------------------|
| **run_noiseq_DE.R**    | Count matrix + sample metadata       | DEGs table                   | For non-parametric differential expression                    | Distribution-free, robust to violations of statistical assumptions |
| **run_bayseq_DE.R**    | Count matrix + sample metadata       | DEGs table                   | Bayesian approach to DE analysis                               | Provides posterior probabilities of differential expression |
| **run_sleuth_DE.R**    | Kallisto quantification output + metadata | DEGs table               | When using pseudoalignment quantifiers (Kallisto/Salmon)       | Works directly on transcript-level abundances, accounts for uncertainty |
| **run_mast_DE.R**      | Expression matrix + metadata         | DEGs table                   | For single-cell RNA-seq, can be shown as advanced option       | Hurdle model handling zero inflation, robust for scRNA-seq data |

---

## 📝 Notes

- **Choice of method depends on data design**:  
  - DESeq2 and edgeR are **default standards**.  
  - limma-voom is ideal for **large datasets**.  
  - NOISeq and baySeq are **alternative statistical frameworks**.  
  - Sleuth connects directly to **pseudoalignment outputs**.  
  - MAST is mostly for **single-cell RNA-seq**, included here as an advanced showcase.  
- **Meta-analysis** is useful when you want a **robust consensus** across multiple methods.  
- Always inspect **diagnostic plots (MA plots, dispersion, PCA, etc.)** before interpreting results.

