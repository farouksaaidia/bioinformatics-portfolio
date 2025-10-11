#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(optparse)
  library(dplyr)
  library(AUCell)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Input Seurat .rds"),
  make_option(c("-m","--marker_sets"), type="character", help="CSV with columns: gene, set_name (marker gene sets)"),
  make_option(c("-o","--output"), type="character", help="Output directory for stats (CSV + RDS)"),
  make_option(c("-k","--top_n"), type="integer", default=50, help="Top N markers per cluster to test (default:50)")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$marker_sets) || is.null(opt$output)) stop("❌ Provide --input --marker_sets --output")

dir.create(opt$output, showWarnings = FALSE, recursive = TRUE)
so <- readRDS(opt$input)
markers <- read.csv(opt$marker_sets, stringsAsFactors = FALSE)

cluster_markers <- FindAllMarkers(so, only.pos = TRUE)
top_markers <- cluster_markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = opt$top_n, with_ties = FALSE)

# Fisher exact: for each cluster and each set, create 2x2 table: in_topN & in_set
results <- data.frame()
for (cl in unique(top_markers$cluster)) {
  genes_top <- top_markers %>% filter(cluster == cl) %>% pull(gene)
  for (set_name in unique(markers$set_name)) {
    set_genes <- markers %>% filter(set_name == set_name) %>% pull(gene)
    a <- sum(genes_top %in% set_genes)        # in top & in set
    b <- length(genes_top) - a               # in top & not in set
    c <- sum(set_genes %in% rownames(so@assays$RNA@counts)) - a  # not in top & in set
    d <- nrow(so) - a - b - c                # not in top & not in set
    tab <- matrix(c(a,c,b,d), nrow=2)
    p <- tryCatch({ fisher.test(tab, alternative = "greater")$p.value }, error=function(e) NA)
    results <- rbind(results, data.frame(cluster=cl, set_name=set_name, overlap=a, p_value=p))
  }
}

write.csv(results, file=file.path(opt$output, "marker_fisher_results.csv"), row.names = FALSE)

# AUCell per set (cell-wise)
expr <- as.matrix(GetAssayData(so, assay="RNA", slot="counts"))
geneSets <- split(markers$gene, markers$set_name)
cellsRankings <- AUCell_buildRankings(expr, nCores=1, plotStats=FALSE)
auc <- AUCell_calcAUC(geneSets, cellsRankings)
auc_df <- as.data.frame(getAUC(auc))
write.csv(auc_df, file=file.path(opt$output, "AUCell_scores.csv"))

saveRDS(list(fisher=results, auc=auc_df), file=file.path(opt$output, "marker_enrichment_stats.rds"))
message("✅ Marker enrichment stats written to ", opt$output)
