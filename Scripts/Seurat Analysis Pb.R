library(tidyverse)
library(Seurat)
library(cowplot)
library(patchwork)
library(magrittr)
library(plotly)
library(viridis)
library(EnhancedVolcano)

setwd("E:/OneDrive - University of Glasgow/2 Projects/2 scRNAseq/S BM B scRNAseq/3 Code")

datadir = "Raw data/"
input.filepath = "Analysis results/1 QC/"
output.filepath = "Analysis results/2 Analysis/"

# Gene list preparation ----
Pb.genes <- readLines(paste0(datadir, "Pb.genes.csv"))

# Loading File ----
  # Requires Seurat object generated by "QC & normalisation & integration"
  
  Pb.combined <- readRDS(paste0(input.filepath, "Pb.combined.rds"))
  head(Pb.combined[[]])
  
  # Host cell binning based on MmUMI ----
  
  Pb.combined$Retic <- cut(Pb.combined$MmUMI, breaks = c(0, 100, 100000), 
                                 labels = c("Normo", "Retic"), 
                                 include.lowest = T)
  
# UMAP and Clustering Parameter optimisation ----
 
  DefaultAssay(Pb.combined) <- "integrated"
  
  DimPlot(Pb.combined, reduction = "pca", group.by = "orig.ident")
  ElbowPlot(Pb.combined, ndims = 50)
  
  for (i in seq(24, 34, 2)) {
  
  PCA.dim <- 1:i
  res <- 0.5
  
  Pb.combined %<>%  
    RunUMAP(reduction = "pca", dims = PCA.dim) %>% 
    FindNeighbors(reduction = "pca", dims = PCA.dim)
  Pb.combined %<>% FindClusters(resolution = res)
  
  print(DimPlot(Pb.combined,reduction = "umap") + ggtitle(paste0("PC 1 to ",i)))
  }
  
  PCA.dim <- 1:24
  
  Pb.combined %<>%  
    RunUMAP(reduction = "pca", dims = PCA.dim) %>% 
    FindNeighbors(reduction = "pca", dims = PCA.dim)
  DimPlot(Pb.combined,reduction = "umap")
  
  for (i in seq(0.2, 1, 0.2)) {
    
    res <- i
    
    Pb.combined %<>% FindClusters(resolution = res)
    
    print(DimPlot(Pb.combined,reduction = "umap") + ggtitle(paste0("res ",i)))
    
  }
  
# 3D UMAP and Clustering ----
   
  PCA.dim <- 1:24
  res <- 0.51
  
  Pb.combined %<>%  
    RunUMAP(reduction = "pca", dims = PCA.dim, n.components = 3) %>% 
    FindNeighbors(reduction = "pca", dims = PCA.dim)
  
  Pb.combined %<>% FindClusters(resolution = res)
  
  DimPlot(Pb.combined,reduction = "umap", dims = c(1,2)) #+ ggtitle(res)
 
   
  p1 <- DimPlot(Pb.combined,reduction = "umap", dims = c(1,2))
  p2 <- DimPlot(Pb.combined,reduction = "umap", dims = c(1,3))
  p3 <- DimPlot(Pb.combined,reduction = "umap", dims = c(2,3))
  p1 + p2 + p3

  Pb.combined$UMAP3D_1 <- Pb.combined[["umap"]]@cell.embeddings[,1]
  Pb.combined$UMAP3D_2 <- Pb.combined[["umap"]]@cell.embeddings[,2]
  Pb.combined$UMAP3D_3 <- Pb.combined[["umap"]]@cell.embeddings[,3]

# 2D UMAP  ----
  
  DefaultAssay(Pb.combined) <- "integrated"
  
  PCA.dim <- 1:24
  res <- 0.51
  
  Pb.combined %<>%  
    RunUMAP(reduction = "pca", dims = PCA.dim) %>% 
    FindNeighbors(reduction = "pca", dims = PCA.dim)
  
  Pb.combined %<>% FindClusters(resolution = res)
  
  DimPlot(Pb.combined, label = TRUE) + NoLegend() 
  
  # Sorting and Renaming stages
  Idents(Pb.combined) <- Pb.combined$seurat_clusters
  order.levels <- c(4,0,1,2,3,5,6,9,7,8,10,11)
  Idents(Pb.combined) <- factor(Idents(Pb.combined), levels = order.levels)
  Pb.combined[["Stages"]] <- Idents(Pb.combined)
  
  heatmap(table(Pb.combined$Stages, Pb.combined$absclust3), Colv = NA, Rowv = NA) 

  Pb.cluster.ids <- c("Early Rings", "Mid Rings", "Late Rings", "Early Trophs 1", 
                      "Early Trophs 2","Mid Trophs", "Late Trophs",
                      "Schizonts", "Outlier", "Early Gams", "Males", "Females")
  names(Pb.cluster.ids) <- levels(Pb.combined)
  Pb.combined <- RenameIdents(Pb.combined, Pb.cluster.ids)
  Pb.combined[["Stages"]] <- Idents(Pb.combined)
  DimPlot(Pb.combined, label = TRUE) + NoLegend()
  DimPlot(subset(Pb.combined, Organ == "MCA"), group.by = "absclust3", label = TRUE) + NoLegend()
  heatmap(table(Pb.combined$Stages, Pb.combined$absclust3), Colv = NA, Rowv = NA) 

  clipr::write_clip(table(Pb.combined$Stages, Pb.combined$orig.ident))
  
# NoMCA subset ----

  Pb.combined.noMCA <- subset(Pb.combined, orig.ident == "SeuratProject", invert = TRUE)
  clipr::write_clip(table(Pb.combined.noMCA$Stages, Pb.combined.noMCA$orig.ident))

# Marker Genes ----

  DefaultAssay(Pb.combined.noMCA) <- "RNA"
  
    # All conserved marker
      ncluster <- nlevels(Pb.combined.noMCA$Stages)
      cluster.list <- vector(mode = "list", length = ncluster)
      names(cluster.list) <- levels(Pb.combined.noMCA$Stages)
      
      for (i in 1:nlevels(Pb.combined.noMCA$Stages)) {
        print(paste0("Working on: #", i, " : ", names(cluster.list)[i], " ..."))
        cluster.list[[i]] <- 
          FindConservedMarkers(Pb.combined.noMCA, ident.1 = names(cluster.list)[i], grouping.var = "orig.ident", 
                               assay = "RNA", verbose = TRUE, min.cells.group = 0) %>% 
          rownames_to_column("GeneID")
        cluster.list[[i]][,"ClusterID"] <- names(cluster.list[i]) 
      }
      
      StagesMarker <- cluster.list

      Pb.combined.marker <- do.call("rbind", StagesMarker[c(1:7, 9:11)])
      Schizonts.marker <- do.call("rbind", StagesMarker[c(8)])
      clipr::write_clip(Pb.combined.marker) 
      clipr::write_clip(Schizonts.marker)
      
# Removal of Outlier cluster ----

  Pb.combined.noMCA <- subset(Pb.combined.noMCA, Stages == "Outlier", invert = T)
  DimPlot(Pb.combined.noMCA)
          
# ADT Analysis ----
      
  DefaultAssay(Pb.combined.noMCA) <- "ADT"
  
  Pb.combined.noMCA$MmUMI1 <- Pb.combined.noMCA$MmUMI + 1

  FeatureScatter(Pb.combined.noMCA, "CD71", "MmUMI1") + scale_y_log10()
  FeatureScatter(Pb.combined.noMCA, "CD44", "MmUMI1") + scale_y_log10()
  
  CD71 <- Pb.combined.noMCA@assays$ADT@data[row.names(Pb.combined.noMCA@assays$ADT@data) %in% "CD71",]
  Pb.combined.noMCA$CD71pos <- cut(CD71, breaks = c(0,1.5,10),
                             labels = c("neg", "pos"), include.lowest = T)
  
  Pb.combined.noMCA$CD71_Retic <- paste0(Pb.combined.noMCA$Retic, "_", Pb.combined.noMCA$CD71pos)
  
  # Percentage cells positive for a marker ----
  
  Pb.Early <- subset(Pb.combined.noMCA, Stages == "Schizonts" | Stages == "Late Trophs" | 
                      Stages == "Mid Trophs" | Stages == "Males" | Stages == "Females", 
                     invert = TRUE)
  Pb.Early.Retic <- subset(Pb.Early, Retic == "Retic")
  Pb.Early.Normo <- subset(Pb.Early, Retic == "Normo")
  
  DefaultAssay(Pb.Early.Retic) <- "RNA"
  DefaultAssay(Pb.Early.Normo) <- "RNA"
  
  # Calculating percent positive per gene, stage and biorep
  
  a <- DotPlot(Pb.Early.Retic, features = unique(Pb.genes), split.by = "orig.ident", cols = viridis(6))
  pct.Retic <- a$data
  pct.Retic$ID <- paste0(pct.Retic$features.plot, "_", pct.Retic$id)

  b <- DotPlot(Pb.Early.Normo, features = unique(Pb.genes), split.by = "orig.ident", cols = viridis(6))
  pct.Normo <- b$data
  pct.Normo$ID <- paste0(pct.Normo$features.plot, "_", pct.Normo$id)
  
  c <- DotPlot(Pb.Early, features = unique(Pb.genes), split.by = "orig.ident", cols = viridis(6))
  pct.all <- c$data
  pct.all$ID <- paste0(pct.all$features.plot, "_", pct.all$id)
 
  pct <- merge(pct.Normo, pct.Retic, by = "ID")
  
  pct$avg.exp.N <- pct$avg.exp.x
  pct$avg.exp.R <- pct$avg.exp.y
  pct$Normo <- pct$pct.exp.x
  pct$Retic <- pct$pct.exp.y
  pct$GeneID <- pct$features.plot.x
  pct$Stage_rep <- pct$id.x
  pct <- pct[,c(14, 15, 16, 17, 18, 19)]
  pct <- separate(pct, Stage_rep, c("Stage", "orig.ident"), sep = "_")
  pct <- separate(pct, orig.ident, c("Organ", "Biorep"), sep = -1)
  
  pct.all$GeneID <- pct.all$features.plot
  pct.all <- separate(pct.all, id, c("Stage", "orig.ident"), sep = "_")
  pct.all <- separate(pct.all, orig.ident, c("Organ", "Biorep"), sep = -1)
  pct.all <- pct.all[,c(1, 2, 3, 4, 5, 6)]
  
# Saving Files ----
  
  saveRDS(Pb.combined, paste0(output.filepath,"Pb.combined.rds"))
  saveRDS(Pb.combined.noMCA, paste0(output.filepath,"Pb.combined.noMCA.rds"))
  
  saveRDS(StagesMarker, paste0(output.filepath,"Pb.ConservedMarker.rds"))

  saveRDS(pct, paste0(output.filepath,"Pb.pct.rds"))
  saveRDS(pct.all, paste0(output.filepath,"Pb.pct.all.rds"))
 