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
Pb.combined.noMCA <- readRDS(paste0(input.filepath, "Pb.combined.NoMCA.rds"))

# Subset to clusters for DEG

Pb.combined.DEG <- subset(Pb.combined.noMCA, Stages == "Schizonts" | Stages == "Late Trophs" |
                            Stages == "Males" | Stages == "Females" | Stages == "Mid Trophs",
                            invert = TRUE)
rm(Pb.combined.noMCA)

# Convert to SingleCellExperiment ----

Pb.sce <- as.SingleCellExperiment(Pb.combined.DEG, assay = "RNA")
plotUMAP(Pb.sce, colour_by = "Stages")
assays(Pb.sce)
dim(counts(Pb.sce))
counts(Pb.sce)[1:6,1:6]
dim(colData(Pb.sce))
head(colData(Pb.sce))

# Simplify Metadata

colData(Pb.sce) %>% as.data.frame %>% 
  transmute(Organ, Biorep, RBC = Retic, ClusterID = Stages, 
            SampleID = paste0(Pb.sce$Organ, Pb.sce$Biorep, Pb.sce$Retic)) %>%
  mutate_all(as.factor) %>% 
  set_rownames(colnames(Pb.sce)) %>% 
  DataFrame -> colData(Pb.sce)
head(colData(Pb.sce))

# Metrics for aggregation ----

# ClusterID (CIDs) and SampleID (SIDs) names and numbers
CIDs <- purrr::set_names(levels(Pb.sce$ClusterID))
nC <- length(CIDs)
SIDs <- purrr::set_names(levels(Pb.sce$SampleID))
nS <- length(SIDs)

# Pseudo-bulk aggregation ----

# Calculate Pseudo-bulk counts by aggregating clusters

groups <- colData(Pb.sce)[,c("ClusterID", "SampleID")]
table(groups)

Pb.agg <- aggregate.Matrix(t(counts(Pb.sce)), groupings = groups, fun = "sum")
Pb.agg[1:6,1:6]

# Vectors identifying main metadata per ClusterID_SampleID
# As not all clusters are present in all bioreps and cell types, this generates
# vectors identifying SampleID, ClusterID, Organ, Biorep and Retic per cell type

splitCID <- stringr::str_split(rownames(Pb.agg), pattern = "_",  n = 2) %>% 
  sapply('[', 1)

splitSID <- stringr::str_split(rownames(Pb.agg), pattern = "_",  n = 2) %>% 
  sapply('[', 2) 

splitRBC <- stringr::str_split(rownames(Pb.agg), pattern = "_",  n = 2) %>% 
  sapply('[', 2) %>% stringr::str_sub(-5, -1)

splitBioRep <- stringr::str_split(rownames(Pb.agg), pattern = "_",  n = 2) %>% 
  sapply('[', 2) %>% stringr::str_sub(-6,-6)

splitOrgan <- stringr::str_split(rownames(Pb.agg), pattern = "_",  n = 2) %>% 
  sapply('[', 2) %>% stringr::str_sub(-8,-7)

# Resplitting matrix by ClusterID and transform

Pb.agg <- split.data.frame(Pb.agg, factor(splitCID)) %>%
  lapply(function(u)  set_colnames(t(u), stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))


# Pseudo-bulk-level MDS Plot ----

mds <- Pb.agg %>% lapply(as.data.frame.matrix) %>% bind_cols %>%
  DGEList(remove.zeros = TRUE) %>%
  calcNormFactors %>%
  plotMDS.DGEList(plot = FALSE)

gg_df <- data.frame(mds[c("x", "y")],
                    ClusterID = splitCID,
                    SampleID = splitSID,
                    Organ = splitOrgan,
                    RBC = splitRBC,
                    Biorep = splitBioRep)

ggplot(gg_df, aes(x, y, col = ClusterID, shape = Biorep)) + 
  geom_point(size = 3, alpha = 0.8) +
  labs(x = "MDS dim. 1", y = "MDS dim. 2") + 
  theme(panel.grid.minor = element_blank()) +
  coord_fixed() + theme_bw()

# DEG using edgeR ----
# Script provided by Dario

pb.agg <- Pb.agg
tt <- list()
allcpm <- list()

for(clst in names(pb.agg)) {
  print(clst)
  dat <- as.matrix(pb.agg[[clst]])
  
  targets <- data.table(
    sample_id = colnames(dat)
  )
  targets[, Organ := ifelse(grepl('^S\\d', sample_id), 'Spleen',
                            ifelse(grepl('^BM\\d', sample_id), 'Bone_marrow',
                                   ifelse(grepl('^B\\d', sample_id), 'Blood', NA)))]
  targets[, Biorep := sub('Normo$|Retic$', '', sub('^B|^BM|^S', '', sample_id))]
  targets[, RBC := ifelse(grepl('Normo$', sample_id), 'normo', 
                          ifelse(grepl('Retic$', sample_id), 'retic', NA))]
  stopifnot(complete.cases(targets))

  design <- model.matrix(~ targets$Biorep + targets$Organ + targets$RBC)
  colnames(design) <- sub('targets$', '', colnames(design), fixed = TRUE)
  colnames(design) <- gsub('\\(|\\)', '', colnames(design))
  
  cnt <- dat[filterByExpr(dat, min.count = 5, min.prop = 0.5), ]
  
  y <- DGEList(cnt)
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  
  fit <- glmQLFit(y, design)
  fit <- glmQLFTest(fit, coef = 5)
  xtt <- as.data.table(topTags(fit, n = Inf, sort.by = "p.value")$table, keep.rownames = 'gene_id')
  xtt[, cluster := clst]
  tt[[length(tt) + 1]] <- xtt

  
  ycpm <- edgeR::cpm(y, log = TRUE, prior.count = 1)
  ycpm <- as.data.table(ycpm, keep.rownames = 'gene_id')
  ycpm <- melt(data = ycpm, id.vars = 'gene_id', variable.name = 'sample_id', value.name = 'logcpm')
  ycpm <- merge(ycpm, targets, by = 'sample_id')
  
  ycnt <- as.data.table(y$counts, keep.rownames = 'gene_id')
  ycnt <- melt(data = ycnt, id.vars = 'gene_id', variable.name = 'sample_id', value.name = 'count')
  ycpm <- merge(ycpm, ycnt, by = c('sample_id', 'gene_id'))   
  ycpm[, cluster := clst]
  allcpm[[length(allcpm) + 1]] <- ycpm
}

tt <- rbindlist(tt)
allcpm <- rbindlist(allcpm)
clipr::write_clip(tt)
clipr::write_clip(allcpm)

# Analyse individual comparisons --- 

clipr::write_clip(tt[, list(n_fdr = sum(FDR < 0.05)), by = list(cluster)][order(cluster)])

# Retic vs Normo ----

Pb.RvsN <- tt
Pb.RvsN <- split(Pb.RvsN, f = Pb.RvsN$cluster)

# Save files

saveRDS(Pb.RvsN, paste0(output.filepath,"Pb.RvsN.rds"))
