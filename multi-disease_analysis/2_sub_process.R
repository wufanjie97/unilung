library(Seurat)
library(tidyverse)
library(cowplot)

dir <- "multi-disease_analysis"


myprocess <- function(adata){
  library(reticulate)
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  library(harmony)
  use_condaenv("/miniconda2/envs/scvi")
  scanpy <- import("scanpy")
  pd <- import("pandas")
  np <- import("numpy")
  
  
  raw_ad <- scanpy$read_h5ad(file.path(dir,paste0(adata,"_raw.h5ad")))
  count <- np$transpose(raw_ad$X)
  colnames(count) <- rownames(raw_ad$obs)
  rownames(count) <- rownames(raw_ad$var)
  mt <- paste0(adata,"_meta.csv")
  meta <- read.csv(file.path(dir,mt))
  rownames(meta) <- colnames(count)
  seu <- CreateSeuratObject(counts = count, meta.data = meta)
  svraw <- paste0(adata,"_raw.rds")
  saveRDS(seu,file.path(dir,svraw))
  print("=============== Done save raw rds ===============")
  
  disease_counts <- table(seu$donor_status)
  ifelse(adata=='mono',cellcut <- 600,ifelse(adata=='bcell',cellcut <- 1400,cellcut <- 2200))
  diseases_to_remove <- names(disease_counts[disease_counts < cellcut])
  diseases_to_keep <- setdiff(seu$donor_status,diseases_to_remove)
  sce <- subset(seu, donor_status %in% diseases_to_keep)
  
  sce <- NormalizeData(sce,normalization.method = "LogNormalize",scale.factor = 10000)
  sce <- FindVariableFeatures(sce,reductionselection.method = "vst",nfeatures = 3000)
  sce <- ScaleData(sce)
  sce <- RunPCA(object = sce,npcs = 30,pc.genes=VariableFeatures(sce),verbose = F)
  sce <- RunHarmony(sce, group.by.vars = "donor_ID", max.iter.harmony = 10)
  sce <- sce %>% 
    RunUMAP(reduction = "harmony", dims = 1:30) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
    FindClusters(resolution = 0.05)
  
  svprocess <- paste0(adata,"_proc.rds")
  saveRDS(sce,file.path(dir,svprocess))
  write.csv(sce@assays$RNA$counts, file.path(dir,paste0(adata,"_proc.csv")), row.names = T)
  
  print("=============== Done UMAP ===============")
  
  svpdf <- paste0(adata,"_celltype.pdf")
  plot_grid(ncol = 3,
            DimPlot(sce,reduction = 'umap',label = F,group.by = "second_anno")+NoAxes(),
            DimPlot(sce,reduction = 'umap',label = F,group.by = "seurat_clusters")+NoAxes(),
            DimPlot(sce,reduction = 'umap',label = F,group.by = "donor_status")+NoAxes())
  ggplot2::ggsave(filename = file.path(dir,svpdf),width = 15,height = 7)
  return(sce)
}


myprocess("mono")
myprocess("bcell")


mono_cell <- readRDS(file.path(dir,'mono_proc.rds'))
mono_cell[["RNA"]]<-as(object=mono_cell[["RNA"]],Class="Assay5")
mono_cell2 <- SketchData(object = mono_cell,ncells = 11000,method = "LeverageScore",sketched.assay = "bpcell")
mono_cell3 <- mono_cell2@assays$bpcell$counts
sce_sub <- subset(mono_cell2,cells=colnames(mono_cell3))
saveRDS(sce_sub,file.path(dir,'mono_sub_raw.rds'))
plot_grid(ncol = 3,
          DimPlot(sce_sub,reduction = 'umap',label = F,group.by = "second_anno")+NoAxes(),
          DimPlot(sce_sub,reduction = 'umap',label = F,group.by = "seurat_clusters")+NoAxes(),
          DimPlot(sce_sub,reduction = 'umap',label = F,group.by = "donor_status")+NoAxes())






bcell <- readRDS(file.path(dir,'bcell_proc.rds'))
bcell[["RNA"]]<-as(object=bcell[["RNA"]],Class="Assay5")
bcell2 <- SketchData(object = bcell,ncells = 12500,method = "LeverageScore",sketched.assay = "bpcell")
bcell3 <- bcell2@assays$bpcell$counts
sce_sub <- subset(bcell2,cells=colnames(bcell3))
saveRDS(sce_sub,file.path(dir,'bcell_sub_raw.rds'))
plot_grid(ncol = 3,
          DimPlot(sce_sub,reduction = 'umap',label = F,group.by = "second_anno")+NoAxes(),
          DimPlot(sce_sub,reduction = 'umap',label = F,group.by = "seurat_clusters")+NoAxes(),
          DimPlot(sce_sub,reduction = 'umap',label = F,group.by = "donor_status")+NoAxes())
