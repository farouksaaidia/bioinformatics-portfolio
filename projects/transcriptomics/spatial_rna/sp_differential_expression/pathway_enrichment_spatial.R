#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
})

option_list <- list(
  make_option(c("-d","--deg_dir"), type="character", help="Directory with DEG CSV files"),
  make_option(c("-o","--output_dir"), type="character", help="Output directory for enrichment CSVs"),
  make_option(c("--fdr_cutoff"), type="double", default=0.05, help="Adjusted p-value cutoff for DEGs (default 0.05)")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$deg_dir) || is.null(opt$output_dir)) stop("❌ Provide --deg_dir and --output_dir")
dir.create(opt$output_dir, recursive=TRUE, showWarnings=FALSE)

files <- list.files(opt$deg_dir, pattern="\\.csv$", full.names=TRUE)
if (length(files) == 0) stop("❌ No DEG CSV files found in directory")

for (f in files) {
  name <- tools::file_path_sans_ext(basename(f))
  degs <- read.csv(f)
  sig <- degs %>% filter(p_val_adj <= opt$fdr_cutoff) %>% pull(gene)
  if (length(sig) == 0) {
    message("⚠️ No significant genes for ", name, "; skipping enrichment")
    next
  }
  # convert to entrez
  entrez <- bitr(sig, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  if (nrow(entrez) == 0) next
  ego <- enrichGO(entrez$ENTREZID, OrgDb = org.Hs.eg.db, ont="BP", readable=TRUE)
  ek <- enrichKEGG(entrez$ENTREZID, organism = "hsa")
  write.csv(as.data.frame(ego), file.path(opt$output_dir, paste0(name, "_GO.csv")), row.names=FALSE)
  write.csv(as.data.frame(ek), file.path(opt$output_dir, paste0(name, "_KEGG.csv")), row.names=FALSE)
  message("✅ Enrichment saved for ", name)
}

message("🎯 Pathway enrichment complete. Results in ", opt$output_dir)
