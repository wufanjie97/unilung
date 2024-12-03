library(Seurat)
library(tidyverse)
library(cowplot)
library(monocle3)
library(harmony)
library(ggsci)
dir <- "lung_cancer_sc/monocle3"
options(Seurat.object.assay.version = "v3")
setwd(dir)

rawsce <- readRDS(file.path(dir,'epi_raw.rds'))
sce <- readRDS(file.path(dir,'epi.rds'))
small <- sce[,sce$annotation=="Malignant cell"]
sce2 <- subset(rawsce,cells=colnames(small))
sce2 <- NormalizeData(sce2,normalization.method = "LogNormalize",scale.factor = 10000)
sce2 <- FindVariableFeatures(sce2,reductionselection.method = "vst",nfeatures = 3000)
sce2 <- ScaleData(sce2)
sce2 <- RunPCA(object = sce2,npcs = 30,pc.genes=VariableFeatures(sce2),verbose = F)
sce2 <- RunHarmony(sce2, group.by.vars = "donor_ID", max.iter.harmony = 10)
sce2 <- sce2 %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = 0.01)

plot_grid(ncol=2,
          DimPlot(sce2,reduction = "umap",cols = pal_lancet("lanonc")(9),
                  label = F,group.by = "seurat_clusters")+ggtitle("Cluster"),
          DimPlot(sce2,reduction = "umap",cols = c("#B24745FF","#79AF97FF","#7B4173FF"),
                  label = F,group.by = "sub_atlas")+ggtitle("Cancer type"))
saveRDS(sce2,file = 'maligant.rds')


sce3 <- SketchData(object = sce2,ncells = 50000,method = "LeverageScore",sketched.assay = "bpcell")
sce4 <- sce3@assays$bpcell@counts
sce_sub <- subset(sce2,cells=colnames(sce4))
sce_sub_tmp <- subset(sce,cells=colnames(sce4))
sce_sub$annotation <- sce_sub_tmp$annotation
sce_sub$stage <- sce_sub_tmp$stage
plot_grid(ncol=2,
          DimPlot(sce_sub,reduction = "umap",cols = pal_lancet("lanonc")(9),
                  label = F,group.by = "seurat_clusters")+ggtitle("Cluster"),
          DimPlot(sce_sub,reduction = "umap",cols = c("#B24745FF","#79AF97FF","#7B4173FF"),
                  label = F,group.by = "sub_atlas")+ggtitle("Cancer type"))
saveRDS(sce_sub,file = 'malignant_sub.rds')




data <- GetAssayData(sce_sub, assay = 'RNA', slot = 'counts')
cell_metadata <- sce_sub@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds,preprocess_method = "PCA",umap.fast_sgd=TRUE,cores = 4)

## load umap embedding from seurat object
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(sce_sub, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed

cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_grid(ncol = 4,
          plot_cells(cds, color_cells_by = "sub_atlas", label_groups_by_cluster=FALSE,label_cell_groups = FALSE,
                     label_leaves=FALSE, label_branch_points=FALSE)+NoAxes()+
            scale_color_manual(values=c(LUAD = "#B24745FF",LUSC = "#79AF97FF",SCLC = "#7B4173FF")),
          plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster=FALSE,label_cell_groups = FALSE, 
                     label_leaves = FALSE,  label_branch_points = FALSE)+NoAxes(),
          plot_cells(cds, color_cells_by = "partition", label_groups_by_cluster=FALSE,label_cell_groups = FALSE, 
                     label_leaves = FALSE,  label_branch_points = FALSE)+NoAxes()+scale_color_npg(),
          plot_cells(cds, color_cells_by = "stage", label_groups_by_cluster=FALSE,label_cell_groups = FALSE, 
                     label_leaves = FALSE,  label_branch_points = FALSE)+NoAxes()+scale_color_jama())
save(sce_sub,cds,file = "malignant_monocle3.RData")




## pseudo gene analysis
load("malignant_monocle3.RData")
cds <- order_cells(cds)

plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster=FALSE,label_cell_groups = FALSE, 
           label_leaves = FALSE,  label_branch_points = FALSE)+NoAxes()

trace_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 1)

track_genes_sig <- trace_genes %>%
  top_n(n=30, morans_I) %>%
  pull(gene_short_name) %>%
  as.character()

choose_gene <- c("TOP2A","MKI67","MUC1","TPX2")
lineage_cds <- cds[rowData(cds)$gene_short_name %in% choose_gene,]

plot_genes_in_pseudotime(lineage_cds,color_cells_by="sub_atlas",min_expr=0.5)+
  scale_color_manual(values=c(LUAD = "#B24745FF",LUSC = "#79AF97FF",SCLC = "#7B4173FF"))
ggplot2::ggsave(filename = "malignant_psd_gene.pdf",width = 6,height = 12)


# umap mapping
plot_cells(cds,genes=choose_gene,
           cell_size=1,
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)
ggplot2::ggsave(filename = "malignant_psd_gene_umap.pdf",width = 8,height = 10)




