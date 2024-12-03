library(Seurat)
library(tidyverse)
library(cowplot)
library(circlize)
library(ComplexHeatmap)

dir <- "multi-disease_analysis"


celltype='mono'
load(paste0(celltype,'_nmf.rdata'))
sce <- readRDS(paste0(celltype,'_sub_raw.rds'))

## nmf heatmap
w = nmf_results@fit@H
colnames(w) = colnames(nmf_counts)
rownames(w) = paste0('NMF', seq(1,ncol(basis(nmf_results))))
w_ = w %>% t %>% scale

col_fun <- colorRamp2(c(-2, 0, 2), c("#8c510a", "white", "#01665e"))
row_ha <- rowAnnotation(donor_status = sce$donor_status,cluster = sce$seurat_clusters)
pdf(paste0(celltype,'_nmf_heatmap.pdf'),width = 8,height = 15)
ComplexHeatmap::Heatmap(w_, 
                        name = "Loading", 
                        border = "black",
                        col = col_fun,
                        width = 4, height = 4, 
                        right_annotation = row_ha,
                        row_km = ncol(basis(nmf_results)), 
                        cluster_rows = T,
                        show_row_names=FALSE,
                        clustering_method_columns = 'single')
dev.off()


nmf_group <- predict(nmf_results) %>% as.data.frame()
nmf_group$group <- paste0('NMF',nmf_group[,1])
nmf_group$sample <- rownames(nmf_group)
sce@meta.data$nmf_group <- nmf_group$group

plot_grid(ncol = 3,
          DimPlot(sce,reduction = 'umap',label = F,group.by = "nmf_group")+NoAxes(),
          DimPlot(sce,reduction = 'umap',label = F,group.by = "seurat_clusters")+NoAxes(),
          DimPlot(sce,reduction = 'umap',label = F,group.by = "donor_status")+NoAxes())


## nmf pattern vlnplot
sce2 <- sce
sce2[["nmf"]] <- CreateAssayObject(counts = w)
DefaultAssay(sce2) <- 'nmf'
color20<-c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', 
           '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
           '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', 
           '#aaffc3', '#808000', '#ffd8b1', '#000075')
plot_grid(ncol = 3,
          VlnPlot(sce2, "NMF1", cols=color20[3:19],group.by = "donor_status", pt.size = 0),
          VlnPlot(sce2, "NMF2", cols=color20[3:19],group.by = "donor_status", pt.size = 0),
          VlnPlot(sce2, "NMF3", cols=color20[3:19],group.by = "donor_status", pt.size = 0))


## get marker of mono-1,-2,-3
# mono-1: NMF2+3----cluster0
# mono-2: NMF1----cluster1+4
# mono-3: NMF2----cluster2+3
new_group <- RenameIdents(sce, "0" = "Mono-1", "1" = "Mono-2",
                          "2" = "Mono-3", "3" = "Mono-3", "4" = "Mono-2"
)
new_group$new_group <- Idents(new_group)
sce$new_group <- Idents(new_group)
saveRDS(sce,'mono_sub.rds')


sce2 <- subset(sce,donor_status != 'NSCLC')
Idents(sce2) <- sce2$new_group


# conversed marker, include monocyte and myeloid marker
VlnPlot(sce2, c("CD14","VCAN","FCN1","LYZ","MARCO"),
        group.by = "new_group", flip = F,stack = T,split.by = "new_group")


# monocyte annotation
select_genes <- c(
  "CD14","S100A8","S100A9",
  "CD79B","CD8A","GZMA","NKG7",
  "FCGR3A","CSF1R","FNDC4","LST1"
)

DotPlot(sce2,features = select_genes,group.by = 'new_group',col.max = 3)+coord_flip()

mono_anno <- RenameIdents(sce2, 
                          "Mono-1" = "CD16 Monocyte", 
                          "Mono-2" = "CD14 Monocyte",
                          "Mono-3" = "Lym-monocyte"
)
mono_anno$mono_anno <- Idents(mono_anno)
sce2$mono_anno <- Idents(mono_anno)
saveRDS(sce2,'mono_sub_latest.rds')

VlnPlot(sce2, features = select_genes,
        group.by = "mono_anno",
        stack = T, 
        sort = T, 
        cols = c("#5E4967","#5F8F6A","#DEA164"),
        split.by =  "mono_anno", 
        flip = TRUE) +
  theme(legend.position = "none")


plot_grid(ncol = 2,
          DimPlot(sce2,reduction = "umap",group.by = "mono_anno") + NoAxes(),
          DimPlot(sce2,reduction = 'umap',group.by = "donor_status",label = F) + NoAxes())



library(UCell)
library(stringr)
library(ggplot2)
library(viridis)
lym_features <- list()
lym_features$Lym <- c('CD79B','CD19','PTPRC','GZMA','CD8A',
                       'CD3D','IL7R','CXCR4','NKG7','CD8B')
marker_score <- AddModuleScore_UCell(sce2, features=lym_features)
a <- colnames(marker_score@meta.data) %>% str_subset("_UCell")
FeaturePlot(marker_score,features = a,order = T, ncol = 1, cols = viridis(256))


library(Seurat)
library(sceasy)
library(reticulate)
sce2[["RNA"]]<-as(object=sce2[["RNA"]],Class="Assay")
sceasy::convertFormat(sce2, from="seurat", to="anndata", assay = "RNA",
                      outFile='mono_sub_latest.h5ad')









celltype='bcell'
load(paste0(celltype,'_nmf.rdata'))
sce <- readRDS(paste0(celltype,'_sub_raw.rds'))

## nmf heatmap
w = nmf_results@fit@H
colnames(w) = colnames(nmf_counts)
rownames(w) = paste0('NMF', seq(1,ncol(basis(nmf_results))))
w_ = w %>% t %>% scale

col_fun <- colorRamp2(c(-2, 0, 2), c("#8c510a", "white", "#01665e"))
row_ha <- rowAnnotation(donor_status = sce$donor_status,cluster = sce$seurat_clusters)
pdf(paste0(celltype,'_nmf_heatmap.pdf'),width = 8,height = 15)
ComplexHeatmap::Heatmap(w_, 
                        name = "Loading", 
                        border = "black",
                        col = col_fun,
                        width = 4, height = 4, 
                        right_annotation = row_ha,
                        row_km = ncol(basis(nmf_results)), 
                        cluster_rows = T,
                        show_row_names=FALSE,
                        clustering_method_columns = 'single')
dev.off()

## nmf-cluster-status umap plot
nmf_group <- predict(nmf_results) %>% as.data.frame()
nmf_group$group <- paste0('NMF',nmf_group[,1])
nmf_group$sample <- rownames(nmf_group)
sce@meta.data$nmf_group <- nmf_group$group

plot_grid(ncol = 3,
          DimPlot(sce,reduction = 'umap',label = F,group.by = "nmf_group")+NoAxes(),
          DimPlot(sce,reduction = 'umap',label = F,group.by = "seurat_clusters")+NoAxes(),
          DimPlot(sce,reduction = 'umap',label = F,group.by = "donor_status")+NoAxes())

## nmf pattern vlnplot
sce2 <- sce
sce2[["nmf"]] <- CreateAssayObject(counts = w)
DefaultAssay(sce2) <- 'nmf'
color20<-c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', 
           '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
           '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', 
           '#aaffc3', '#808000', '#ffd8b1', '#000075')
plot_grid(ncol = 4,
          VlnPlot(sce2, "NMF1", cols=color20[6:19],group.by = "donor_status", pt.size = 0),
          VlnPlot(sce2, "NMF2", cols=color20[6:19],group.by = "donor_status", pt.size = 0),
          VlnPlot(sce2, "NMF3", cols=color20[6:19],group.by = "donor_status", pt.size = 0),
          VlnPlot(sce2, "NMF4", cols=color20[6:19],group.by = "donor_status", pt.size = 0))


## get marker of bcell-1,-2,-3
# bcell-1: NMF1+3----cluster0+3+4
# bcell-2: NMF4----cluster1
# bcell-3: NMF1+2+3+4----cluster2
new_group <- RenameIdents(sce, "0" = "Bcell-1", "1" = "Bcell-2",
                          "2" = "Bcell-3", "3" = "Bcell-1", "4" = "Bcell-1"
)
new_group$new_group <- Idents(new_group)
sce$new_group <- Idents(new_group)
saveRDS(sce,'bcell_sub.rds')

## identification of monocyte subset
sce2 <- subset(sce,donor_status != 'NSCLC')
Idents(sce2) <- sce2$new_group


# conversed marker, include monocyte and myeloid marker
VlnPlot(sce2, c("PTPRC","CCL5","CD3D","CD80","CD180","CD37","CD79B","MS4A1","CD74"),
        group.by = "new_group", flip = F,stack = T,split.by = "new_group")


# monocyte annotation
select_genes <- c(
  "CD79B","CD80","CD180",
  "CD3D","CD3E","IL7R","CD3G",
  "CD79A","MS4A1","CD83"
)
DotPlot(sce2,features = unique(select_genes),group.by = 'new_group',col.max = 3)+coord_flip()

bcell_anno <- RenameIdents(sce2, 
                          "Bcell-1" = "Memory B cell", 
                          "Bcell-2" = "T-like B cell",
                          "Bcell-3" = "Activated B cell"
)
bcell_anno$bcell_anno <- Idents(bcell_anno)
sce2$bcell_anno <- Idents(bcell_anno)
saveRDS(sce2,'bcell_sub_latest.rds')


VlnPlot(sce2, features = select_genes,
        group.by = "bcell_anno",
        stack = T, 
        sort = T, 
        cols = c("#5E4967","#5F8F6A","#DEA164"),
        split.by =  "bcell_anno", 
        flip = TRUE) +
  theme(legend.position = "none")


plot_grid(ncol = 2,
          DimPlot(sce2,reduction = "umap",group.by = "bcell_anno") + NoAxes(),
          DimPlot(sce2,reduction = 'umap',group.by = "donor_status",label = F) + NoAxes())



library(UCell)
library(stringr)
library(ggplot2)
library(viridis)
tcell_features <- list()
tcell_features$tcell <- c('CD3E','CD3G','CD8A','CXCL13','GATA3',
                      'CD3D','IL7R','IL2RB','CD8B','CD4')
marker_score <- AddModuleScore_UCell(sce2, features=tcell_features)
a <- colnames(marker_score@meta.data) %>% str_subset("_UCell")
FeaturePlot(marker_score,features = a,order = T, ncol = 1, cols = viridis(256))



library(Seurat)
library(sceasy)
library(reticulate)
sce2[["RNA"]]<-as(object=sce2[["RNA"]],Class="Assay")
sceasy::convertFormat(sce2, from="seurat", to="anndata", assay = "RNA",
                      outFile='bcell_sub_latest.h5ad')