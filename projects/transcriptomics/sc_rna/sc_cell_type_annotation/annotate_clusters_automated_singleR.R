#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(SingleR)
  library(optparse)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Input Seurat .rds file (raw or clustered)"),
  make_option(c("-o", "--output"), type="character", help="Output annotated Seurat .rds file"),
  make_option(c("-r", "--reference"), type="character", default="HPA", help="Reference to use: HPA | BlueprintEncode | Monaco | custom (path to SCE). Default: HPA"),
  make_option(c("-k", "--assay"), type="character", default="RNA", help="Assay to use for SingleR (default: RNA)")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$output)) stop("❌ Provide both --input and --output")

message("📂 Loading Seurat object...")
seurat_obj <- readRDS(opt$input)

# prepare SingleCellExperiment
sce <- tryCatch({
  as.SingleCellExperiment(seurat_obj, assay = opt$assay)
}, error = function(e) stop("❌ Failed to convert Seurat to SingleCellExperiment: ", e$message))

message("📦 Loading reference...")
ref_choice <- tolower(opt$reference)
if (ref_choice %in% c("hpa", "humanprimarycelldata", "humanprimarycelldata")) {
  if (!requireNamespace("celldex", quietly = TRUE)) {
    stop("❌ celldex required for HPA reference. Install with BiocManager::install('celldex')")
  }
  ref <- celldex::HumanPrimaryCellAtlasData()
} else if (ref_choice %in% c("blueprintencode", "blueprint")) {
  ref <- celldex::BlueprintEncodeData()
} else if (ref_choice %in% c("monaco", "monacoimmune")) {
  ref <- celldex::MonacoImmuneData()
} else if (file.exists(opt$reference)) {
  ref <- readRDS(opt$reference)
} else {
  stop("❌ Unknown reference. Provide 'HPA', 'BlueprintEncode', 'Monaco' or a path to a saved SingleCellExperiment reference.")
}

message("🧠 Running SingleR...")
pred <- SingleR(test = sce, ref = ref, labels = ref$label.main %||% colData(ref)[,1])
# attach labels
pred_labels <- pred$labels
seurat_obj$predicted_cell_type_singleR <- pred_labels[ colnames(sce) ]

saveRDS(seurat_obj, opt$output)
message("✅ SingleR annotation saved to ", opt$output)
