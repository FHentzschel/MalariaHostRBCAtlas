library(tidyverse)
library(magrittr)
library(Matrix.utils)
library(data.table)

library(Seurat)
library(edgeR)
library(scater)

library(patchwork)
library(cowplot)
library(EnhancedVolcano)
library(UpSetR)

setwd("E:/OneDrive - University of Glasgow/2 Projects/2 scRNAseq/S BM B scRNAseq/3 Code")

datadir = "Raw data/"
input.filepath = "Analysis results/2 Analysis/"
output.filepath = "Analysis results/3 DEG/"

# I roughly followed the workflows delineated in these vignettes:
# http://biocworkshops2019.bioconductor.org.s3-website-us-east-1.amazonaws.com/page/muscWorkshop__vignette/
# https://hbctraining.github.io/scRNA-seq/lessons/pseudobulk_DESeq2_scrnaseq.html
# The idea is to do pseudobulk aggregation to be able to apply edgeR for differential gene expression analysis
# Data preparation ----

# Loading file ----
BM.combined <- readRDS(paste0(input.filepath, "BM.combined.rds"))


# Convert to SingleCellExperiment ----

sce <- as.SingleCellExperiment(BM.combined, assay = "RNA")
plotUMAP(sce, colour_by = "BMIdent")
assays(sce)
dim(counts(sce))
counts(sce)[1:6,1:6]
dim(colData(sce))
head(colData(sce))

# Simplify Metadata

colData(sce) %>% as.data.frame %>% 
  transmute(Biorep, Status, ClusterID = BMIdent, 
            SampleID = paste0(sce$Biorep, sce$Status)) %>%
  mutate_all(as.factor) %>% 
  set_rownames(colnames(sce)) %>% 
  DataFrame -> colData(sce)
head(colData(sce))

# Metrics for aggregation ----

# ClusterID (CIDs) and SampleID (SIDs) names and numbers
CIDs <- purrr::set_names(levels(sce$ClusterID))
nC <- length(CIDs)
SIDs <- purrr::set_names(levels(sce$SampleID))
nS <- length(SIDs)

# Pseudo-bulk aggregation ----

# Calculate Pseudo-bulk counts by aggregating clusters

groups <- colData(sce)[,c("ClusterID", "SampleID")]
table(groups)

BM.agg <- aggregate.Matrix(t(counts(sce)), groupings = groups, fun = "sum")
BM.agg[1:6,1:6]

# Resplitting matrix by ClusterID and transform

BM.agg <- split.data.frame(BM.agg, factor(CIDs)) %>%
  lapply(function(u)  set_colnames(t(u), SIDs))

# DEG using edgeR ----
# Script provided by Dario

tt <- list()

for(clst in names(BM.agg)) {
  print(clst)
  dat <- as.matrix(BM.agg[[clst]])
  
  targets <- data.table(
    sample_id = colnames(dat)
  )

  targets[, Biorep := sub('Infected$|Naive$', '', sample_id)]
  targets[, Status := ifelse(grepl('Infected$', sample_id), 'Infected', 
                          ifelse(grepl('Naive$', sample_id), 'Naive', NA))]
  stopifnot(complete.cases(targets))

  design <- model.matrix(~ 0 + targets$Status)
  colnames(design) <- sub('targets$', '', colnames(design), fixed = TRUE)
  colnames(design) <- gsub('\\(|\\)', '', colnames(design))
  
  cnt <- dat[filterByExpr(dat, min.count = 5, min.prop = 0.5), ]
  
  y <- DGEList(cnt)
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)

  fit <- glmQLFit(y, design)
  fit <- glmQLFTest(fit, contrast = c(1,-1))
  xtt <- as.data.table(topTags(fit, n = Inf, sort.by = "p.value")$table, keep.rownames = 'gene_id')
  xtt[, cluster := clst]
  tt[[length(tt) + 1]] <- xtt

}

tt <- rbindlist(tt)
clipr::write_clip(tt)

# Analyse individual comparisons --- 

clipr::write_clip(tt[, list(n_fdr = sum(FDR < 0.05)), by = list(cluster)][order(cluster)])
clipr::write_clip(tt[, list(n_fdr = sum(FDR < 0.1)), by = list(cluster)][order(cluster)])

# Retic vs Normo ----

BM.InfvsNaive <- split(tt, f = tt$cluster)

# Save files

saveRDS(BM.InfvsNaive, paste0(output.filepath,"BM.InfvsNaive.rds"))
