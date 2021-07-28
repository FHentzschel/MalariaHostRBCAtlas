library(dplyr)
library(tibble)
library(Seurat)
library(cowplot)
library(patchwork)
library(scran)
library(scater)
library(magrittr)

setwd("E:/OneDrive - University of Glasgow/2 Projects/2 scRNAseq/S BM B scRNAseq/3 Code")

datadir = "Raw data/"
analysis.filepath = "Analysis results/1 QC/"
referencedir = "Raw data/Reference Datasets/"

flag_calc_visualisation = TRUE

# 1) *** Data preparation *** ----
  # Gene list preparation ----
    Pb.genes <- readLines(paste0(datadir, "Pb.genes.csv"))
    Mm.genes <- readLines(paste0(datadir, "Mm.genes.csv"))
   
  # Loading Data into a list of seurat objects----
    config <- data.frame(
      sample.name = c("S1","S2", "S3", "BM1", "BM2", "BM3", "B1", "B2", "B3")
      ,data.subdir = c("S1","S2", "S3", "BM1", "BM2", "BM3", "B1", "B2", "B3") 
      ,organ.name = c("S", "S", "S", "BM", "BM", "BM", "B", "B", "B")
      ,status.name = c("Infected", "Infected", "Naive", "Infected", "Infected", "Naive", "Infected", "Infected", "Naive")
      ,biorep = c("1", "2", "3", "1", "2", "3", "1", "2", "3")
    )
    
    config$sample.name <- as.character(config$sample.name)
    
    All.list = list()
    for (i in 1:nrow(config)) {
      print(paste("loading dataset: ", config$sample.name[i]))
      counts <- Read10X(data.dir = paste0(datadir, config$data.subdir[i], "/raw_feature_bc_matrix short/"))
      seurat.object <- CreateSeuratObject(counts = counts, project = config$sample.name[i], min.cells = 3, min.features = 100)
      seurat.object <- subset(seurat.object, subset = nCount_RNA > 99)
      print(dim(seurat.object))
      seurat.object[["Organ"]] <- config$organ.name[i]
      seurat.object[["Status"]] <- config$status.name[i]
      seurat.object[["Biorep"]] <- config$biorep[i] 
      All.list[[config$sample.name[i]]] <- seurat.object
    }
    
    rm(seurat.object, counts)

  # Metadata addition ----
    for (listname in names(All.list)) {
      So = All.list[[listname]]  
     
      So[["PbPercent.mt"]] <- PercentageFeatureSet(So, pattern = "*PBANKA-MIT")
      So[["MmPercent.mt"]] <- PercentageFeatureSet(So, pattern = "mt-")
      
      So[["PbUMI"]] <- PercentageFeatureSet(So, pattern = "PBANKA") * So$nCount_RNA/100
      So[["Percent.PbUMI"]] <- PercentageFeatureSet(So, pattern = "PBANKA")
      So[["MmUMI"]] <- So$nCount_RNA - So$PbUMI
     
      So[["PbGenes"]] <- Matrix::colSums(So@assays$RNA@data[row.names(So@assays$RNA@data) %in% Pb.genes,] != 0)
      So[["MmGenes"]] <- So$nFeature_RNA - So$PbGenes
    
      All.list[[listname]] <- So
    }
    rm(So)
    
  # Plotting features, total cell numbers ----
  if (flag_calc_visualisation == TRUE) {
    for (listname in names(All.list)) {
      So = All.list[[listname]]  
      
      p1 <- FeatureScatter(So, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.1)
      p2 <- FeatureScatter(So, feature1 = "PbPercent.mt", feature2 = "MmPercent.mt", pt.size = 0.1)
      p3 <- FeatureScatter(So, feature1 = "PbGenes", feature2 = "MmGenes", pt.size = 0.1)
      p4 <- FeatureScatter(So, feature1 = "PbUMI", feature2 = "MmUMI", pt.size = 0.1)
      p5 <- FeatureScatter(So, feature1 = "MmUMI", feature2 = "MmGenes", pt.size = 0.1)
      p6 <- FeatureScatter(So, feature1 = "PbUMI", feature2 = "PbGenes", pt.size = 0.1)
      
      print(p1 + p2 + p3 + p4 + p5 + p6)
      
      print(listname)
      print(dim(So))
    }
    rm(So)
    
  } # end if flag_calc_visualisation
    

  # Filtering Seurat objects ----
  MmGenes_max = c(3000, 4000, 4000, 7000, 6000, 5000, 200, 100, 250)
  MmUMI_max =   c(20000, 20000, 20000, 50000, 25000, 28000, 14000, 2000, 15000)
  PbGenes_max = c(3000, 3100, 10, 3000, 3100, 10, 3000, 2500, 100)
  PbUMI_max =   c(9000, 20000, 10, 10000, 18000, 10, 12000, 8000, 100)
  
  # Cell numbers before filtering
  
   for (i in seq_along(All.list)) {
     listname = names(All.list)[i]
     
    print(paste("Cell numbers before filtering for dataset: ", listname))
    print(dim(All.list[[listname]]))
    
    All.list[[listname]] %<>% 
      subset(subset = nFeature_RNA > 100 & MmGenes < MmGenes_max[[i]] & MmUMI < MmUMI_max[[i]] &
             PbGenes < PbGenes_max[[i]] & PbUMI < PbUMI_max[[i]] & 
            PbPercent.mt < 1 & MmPercent.mt < 5)
    # Cell numbers after filtering
    print(paste("Cell numbers after filtering for dataset: ", listname))
    print(dim(All.list[[listname]]))
  }
  
  # Adding CITE-Seq Data ----
  
  names.ADT = names(All.list)

  for (i in which(names(All.list) %in% names.ADT)) {
    listname = config$sample.name[i]
    print(listname)
    Pb = All.list[[listname]]
    cs.counts <- Read10X(data.dir = paste0(datadir,"CiteSeq/CiteSeq.RAWBAR.", config$data.subdir[i], "/umi_count short/"), gene.column = 1)
    cs.counts <- cs.counts[setdiff(rownames(cs.counts), c("unmapped")),
                           colnames(cs.counts) %in% colnames(Pb)]
    CS_colnames <- cs.counts %>%
      CreateSeuratObject(project = "CS", min.cells = 3, min.features = 0) %>%
      colnames()

    print(paste("Dims *before* cite seq addition for sample ", listname, ":", nrow(Pb), ncol(Pb)))
    Pb %<>%  subset(cells = CS_colnames)
    Pb[["ADT"]] <- CreateAssayObject(counts = cs.counts)
    print(paste("Dims *after* cite seq addition for sample ", listname, ":", nrow(Pb), ncol(Pb)))
    All.list[[listname]] <- Pb
  }
  rm(cs.counts, Pb, CS_colnames)
  ADT.genes <- rownames(All.list[["B1"]][["ADT"]])
  
  # ADT normalisation and scaling ----
  
    for (i in which(names(All.list) %in% names.ADT)) {
      listname = config$sample.name[i]
      print(listname)
      Pb = All.list[[listname]]
      Pb %<>%  
        NormalizeData(assay = "ADT", normalization.method = "CLR") %>% 
        ScaleData(assay = "ADT")
      All.list[[listname]] <- Pb
    }
  
  # Renaming cell names for unique names----
  
    for (i in 1:length(All.list)) {
      All.list[[i]] <- RenameCells(All.list[[i]], new.names = paste0(colnames(All.list[[i]]), "_", i))
    }
  
  # Saving file ----
    saveRDS(All.list, file = paste0(analysis.filepath, "All.list.rds"))

# 2) *** Mm Analysis *** ----   
  # Mm Dataset selection and subsetting----
    Mm.Datasets <- c("S1", "S2", "S3", "BM1", "BM2", "BM3")
    Mm.list <- All.list[Mm.Datasets]
    for (listname in names(Mm.list)) {
      So <- Mm.list[[listname]]
      print(listname)
      So %>% dim() %>% print()
      So %<>%
        subset(subset = (MmGenes > 100)) %>%
        subset(features = c(Mm.genes, ADT.genes)) 
      So %>% dim() %>% 
        print()
      Mm.list[[listname]] <- So
    }
    rm(So)
    
  # Scran normalisation ----
  
    scran.cluster <- c(100, 300, 100, 100, 100, 100)
    scran.method <- c("hclust", "igraph", "hclust", "hclust", "hclust", "hclust")
    names(scran.cluster) <- names(Mm.list)
    names(scran.method) <- names(Mm.list)
    
    for (listname in names(Mm.list)) {
      So <- Mm.list[[listname]]
      
     ADT.counts <- So %>% GetAssay("ADT")
    
      So %<>% as.SingleCellExperiment()
      qclust <- So %>% scran::quickCluster(min.size = scran.cluster[[listname]], method = scran.method[[listname]])
      So %<>% scran::computeSumFactors(clusters = qclust)
      So %>% sizeFactors() %>% summary() %>% 
        print()
      So %<>% 
        scater::logNormCounts(log = FALSE)
      logcounts(So) <- as.sparse(log(assay(So, "normcounts") + 1))
      So <- as.Seurat(So, counts = "counts", data = "logcounts")
      
     So@assays$ADT <- ADT.counts
  
      Mm.list[[listname]] <- So
    }
    
    rm(So)
    rm(ADT.counts)
    
  # Identification of variable genes ----
  
    for (listname in names(Mm.list)) {
      So <- Mm.list[[listname]]
      So %<>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
      So %>% VariableFeaturePlot(pt.size = 0.1) %>% 
        print()
      Mm.list[[listname]] <- So
    }
    rm(So)
    
  # Saving file ----
    saveRDS(Mm.list, file = paste0(analysis.filepath, "Mm.list.rds")) 
    
  # Performing CCA and merging data sets ----
    Mm.anchors <- FindIntegrationAnchors(Mm.list) 
    Mm.combined <- IntegrateData(Mm.anchors)
    
    rm(Mm.anchors)
    
    DefaultAssay(Mm.combined) <- "integrated"
  
  # Run the standard workflow for visualization and clustering ----
    Mm.combined %<>%
      ScaleData(do.scale = TRUE, features = Mm.genes) %>% 
      RunPCA(npcs = 50, verbose = TRUE)
    DimPlot(Mm.combined, reduction = "pca", group.by = "orig.ident")
    ElbowPlot(Mm.combined, ndims = 50)
    
  # Transfer labels from Bacchin 2020 ----
    # Load Baccin files 
      load(paste0(referencedir, "/Bone Marrow Bacchin 2020/NicheData10x.rda"))
      NicheData10x[["Ident"]] <- Idents(NicheData10x)
  
    # Transfer labels 
      DefaultAssay(Mm.combined) <- "RNA"
      
      BM.anchors <- FindTransferAnchors(reference = NicheData10x, 
                                        query = Mm.combined, normalization.method = "LogNormalize",
                                        dims = 1:30, query.assay = "RNA", reference.assay = "RNA") 
      BM.predictions <- TransferData(anchorset = BM.anchors, refdata = NicheData10x$Ident, dims = 1:30)
      Mm.combined %<>% AddMetaData(metadata = BM.predictions)
      
      VlnPlot(Mm.combined, "prediction.score.max", group.by = "predicted.id", pt.size = 0.1)
    
    # Remove unnecessary files
      rm("NicheData10x", "BM.anchors", "BM.predictions")
      
  # Saving file ----
    saveRDS(Mm.combined, file = paste0(analysis.filepath, "Mm.combined.rds"))

  # Subsetting to original identities ----
  
   for (listname in names(Mm.list)) {
      Mm.list[[listname]] <- Mm.combined %>% subset(subset = (orig.ident == listname))
      print(paste0("Subsetting ", listname, " with orig.ident == ", "'", listname, "'"))
   }
      saveRDS(Mm.list, file = paste0(analysis.filepath, "Mm.list.rds"))
      
      rm(Mm.combined)
  
  # Transfer Mm cluster to main Seurat ----
    # This will add the cluster IDs of the host cell clusters identified in the above Seurat object to the metadata of
    # the main Seurat object
    # Replace NA for the cells not clustered in PbMm by "RBCs"
      
      for (listname in names(All.list)) {
        So = All.list[[listname]]
        
        if (listname %in% names(Mm.list)) {
          Mm = Mm.list[[listname]]
    
           So[["MmIdent"]] <- data.frame(Mm$predicted.id)
           levels(So$MmIdent)
           So$MmIdent <- factor(So$MmIdent, levels = c(levels(So$MmIdent), "RBCs"))
           So$MmIdent[is.na(So$MmIdent)] = "RBCs"
        } else {
           So[["MmIdent"]] <- factor("RBCs")
         }
         print(listname)
         print(levels(So$MmIdent))
         
         All.list[[listname]] <- So
      }
      rm(So)
      rm(Mm)
      rm(Mm.list)
  
  # Saving file ----
    saveRDS(All.list, file = paste0(analysis.filepath, "All.list.rds"))
      
# 3) *** PbAnalysis *** ----
  # Pb Dataset selection and subsetting ----
  
    Pb.Datasets <- c("S1","S2", "BM1","BM2", "B1","B2")
    Pb.list <- All.list[Pb.Datasets]
    rm(All.list)
  
    for (listname in names(Pb.list)) {
      So <- Pb.list[[listname]]
      So %>% dim() %>% 
        print()
      So %<>% 
        subset((MmIdent == "RBCs") | (MmIdent == "Erythroblasts")) %>%
        subset(subset = (PbGenes > 100)) %>% 
        subset(features = c(Pb.genes, ADT.genes))
      So %>% dim() %>% 
        print()
      Pb.list[[listname]] <- So
    }
    rm(So)
    
    # Addition of MCA to Pb.list ----
      # Prepare MCA files 
        MCA.metadata <- read.table(paste0(referencedir,"Malaria Cell Atlas/10X Data/allpb10x/allpb10x_pheno.csv"),
                                   sep = ",", header = TRUE, row.names = 1) %>% 
          as.data.frame()
        
        MCA <- read.table(paste0(referencedir,"Malaria Cell Atlas/10X Data/allpb10x/allpb10x_counts .1.txt"),
                          sep = ",", header = TRUE, row.names = 1) %>% 
          as.matrix() %>% 
          CreateSeuratObject(meta.data = MCA.metadata, project = "MCA", min.cells = 3, min.features = 0)
        
        rm(MCA.metadata)
        
        MCA[["PbPercent.mt"]] <- PercentageFeatureSet(MCA, pattern = "*PBANKA-MIT*")
        MCA[["PbUMI"]] <- MCA$nCount_RNA
        MCA[["PbGenes"]] <- MCA$nFeature_RNA
        MCA[["Organ"]] <- "MCA"
    
      # Subsetting features
        MCA %<>%  subset(features = Pb.genes)
    
      # Addition MCA to Pb.list
        Pb.list.names <- names(Pb.list)
        Pb.list <- c(Pb.list, MCA)
        names(Pb.list) <- c(Pb.list.names,"MCA")
        rm(MCA)
    
    # Scran normalisation ----

      for (listname in names(Pb.list)) {
        So <- Pb.list[[listname]]
        
        if (listname %in% names.ADT) {
         ADT.counts <- GetAssay(So, "ADT")
        }
        
        So %<>% as.SingleCellExperiment()
        qclust <- So %>% scran::quickCluster(min.size = 100, method = "hclust")
        So %<>% scran::computeSumFactors(clusters = qclust)
        So %>% sizeFactors() %>% summary() %>% 
          print()
        So %<>% 
          scater::logNormCounts(log = FALSE)
        logcounts(So) <- as.sparse(log(assay(So, "normcounts") + 1))
        So <- as.Seurat(So, counts = "counts", data = "logcounts")
        
        if (listname %in% names.ADT) {
          So@assays$ADT <- ADT.counts
        }
        
        Pb.list[[listname]] <- So
      }
      rm(So)
      rm(ADT.counts)
    
    # Identification of variable genes ----
      for (listname in names(Pb.list)) {
        So <- Pb.list[[listname]]
        So %<>% 
          FindVariableFeatures(selection.method = "vst", nfeatures = 1000)
        So %>% VariableFeaturePlot(pt.size = 0.1) %>% 
          print()
        Pb.list[[listname]] <- So
      }
      rm(So)
      
  # Saving file ----
      saveRDS(Pb.list, file = paste0(analysis.filepath, "Pb.list.rds"))

  # Performing integration ----
     Pb.combined <- Pb.list %>% 
        FindIntegrationAnchors(normalization.method = "LogNormalize") %>% 
        IntegrateData(normalization.method = "LogNormalize")
     rm(Pb.list)
    
    DefaultAssay(Pb.combined) <- "integrated"
    
  # Scale Data and run PCA ----
    Pb.combined %<>%  
      ScaleData(do.scale = TRUE, features = Pb.genes) %>% 
      RunPCA(npcs = 50, verbose = TRUE)
    
     DimPlot(Pb.combined, reduction = "pca")
    ElbowPlot(Pb.combined, ndims = 50)
    
    # Saving file ----
    saveRDS(Pb.combined, file = paste0(analysis.filepath, "Pb.combined.rds"))
