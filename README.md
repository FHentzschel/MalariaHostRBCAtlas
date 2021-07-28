# Malaria Host RBC Atlas

Raw count matrices and R scripts to analyse scRNA-seq data of <I>Plasmodium berghei</I> parasite cells and murine host cells isolated from spleen, blood and bone marrow. 

# Script guidance

Scripts are to be executed in the following order

1) QC & Normalisation & Integration   
   Generates three integrated seurat objects for  <I>P. berghei</I> cells, spleen cells, and host cells. The following reference data sets were used for annotation of clusters: The malaria cell atlas (Howick, V. M. et al., Science 2019, DOI: 10.1126/science.aaw2619) for  <I>P. berghei</I> cells and the bone marrow niche atlas (Baccin, C. et al., Nature Cell Biology 2020, DOI: 10.1038/s41556-019-0439-6) for host cells.
   
2) Seurat Analysis Pb/BM/S  
   Takes the merged/integrated seurat objects as input and performs UMAP and clustering, as well as identification of marker genes for each cluster. Separate scripts for P.          berghei, spleen and bone marrow 
   
3) DEG edgeR  
   Performs differential gene expression (DGE) analysis using edgeR on the Seurat objects generated in (2). Scripts S and BM calculate DGE between infected and uninfected samples in spleen and bone marrow, respectively, Pb Organs calculates DGE in <I>P. berghei</I> cells isolated from different organs, and Pb Retic calculates DGE in  <I>P. berghei</I> cells between normocyte and reticulocyte host cells. 
   
4) Plots  
   Code to generate plots that are depicted in the manuscript. 
   
