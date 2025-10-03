#!/usr/bin/env Rscript

# Compare DESeq2, edgeR, limma in one run
suppressMessages(library(DESeq2))
suppressMessages(library(edgeR))
suppressMessages(library(limma))

args <- commandArgs(trailingOnly = TRUE)
counts <- read.table(args[1], header=TRUE, row.names=1)
colData <- read.table(args[2], header=TRUE, row.names=1)

# DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ condition)
dds <- DESeq(dds)
res_deseq2 <- results(dds)

# edgeR
group <- factor(colData$condition)
y <- DGEList(counts=counts, group=group)
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit)
res_edger <- topTags(qlf, n=Inf)

# limma-voom
y <- DGEList(counts=counts, group=group)
y <- calcNormFactors(y)
design <- model.matrix(~group)
v <- voom(y, design, plot=FALSE)
fit <- lmFit(v, design)
fit <- eBayes(fit)
res_limma <- topTable(fit, number=Inf, sort.by="P")

write.csv(as.data.frame(res_deseq2), paste0(args[3], "_DESeq2.csv"))
write.csv(as.data.frame(res_edger), paste0(args[3], "_edgeR.csv"))
write.csv(as.data.frame(res_limma), paste0(args[3], "_limma.csv"))
