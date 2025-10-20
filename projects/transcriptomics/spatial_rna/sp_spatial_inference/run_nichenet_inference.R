#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(nichenetr)
  library(dplyr)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Input annotated Seurat .rds"),
  make_option(c("-s","--sender"), type="character", help="Sender cell type(s), comma-separated"),
  make_option(c("-r","--receiver"), type="character", help="Receiver cell type"),
  make_option(c("-o","--output"), type="character", help="Output CSV for ligand activity scores")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$sender) || is.null(opt$receiver) || is.null(opt$output)) {
  stop("❌ Provide input, sender, receiver, and output paths")
}

seurat_obj <- readRDS(opt$input)
sender_ids <- unlist(strsplit(opt$sender, ","))
receiver_id <- opt$receiver

cat("🧩 Running NicheNet for sender(s):", paste(sender_ids, collapse=", "), "-> receiver:", receiver_id, "\n")

# Prepare data
Idents(seurat_obj) <- "predicted_cell_type"
exprs_sender <- AverageExpression(seurat_obj, group.by="predicted_cell_type")$RNA[, sender_ids, drop=FALSE]
exprs_receiver <- AverageExpression(seurat_obj, group.by="predicted_cell_type")$RNA[, receiver_id, drop=FALSE]

# Load ligand-target matrix
ligand_target_matrix <- readRDS(system.file("extdata", "ligand_target_matrix.rds", package="nichenetr"))
lr_network <- readRDS(system.file("extdata", "lr_network.rds", package="nichenetr"))

# Identify potential ligands expressed in senders
expressed_ligands <- rownames(exprs_sender)[rowMeans(exprs_sender) > 0.1]
expressed_receptors <- rownames(exprs_receiver)[rowMeans(exprs_receiver) > 0.1]

# Run ligand activity
ligand_activities <- predict_ligand_activities(geneset=expressed_receptors,
                                               background_expressed_genes=rownames(exprs_receiver),
                                               ligand_target_matrix=ligand_target_matrix,
                                               potential_ligands=expressed_ligands)

write.csv(ligand_activities, opt$output, row.names=FALSE)
cat("✅ NicheNet ligand–receptor inference complete ->", opt$output, "\n")
