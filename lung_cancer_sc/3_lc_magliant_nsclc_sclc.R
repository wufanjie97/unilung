library(Seurat)
library(tidyverse)
library(cowplot)
library(ggsci)
library(harmony)
library(clutree)
dir <- "lung_cancer_sc/nsclc_sclc"
options(Seurat.object.assay.version = "v3")
setwd(dir)

load("malignant_monocle3.RData")# sce=sce_sub

rawsce <- readRDS('epi_raw.rds')
sce <- readRDS('malignant_sub.rds')
subcds <- cds[,cds@clusters$UMAP$partitions %in% c(1,2)]
nsclc_sclc <- rawsce[,colnames(subcds)]
nsclc_sclc_tmp <- sce[,colnames(subcds)]
nsclc_sclc$stage <- nsclc_sclc_tmp$stage
saveRDS(nsclc_sclc,'nsclc_sclc_raw.rds')


sce <- nsclc_sclc
sce2 <- NormalizeData(sce,normalization.method = "LogNormalize",scale.factor = 10000)
sce2 <- FindVariableFeatures(sce2,reductionselection.method = "vst",nfeatures = 3000)
sce2 <- ScaleData(sce2)
sce2 <- RunPCA(object = sce2,npcs = 30,pc.genes=VariableFeatures(sce2),verbose = F)
sce2 <- RunHarmony(sce2, group.by.vars = "donor_ID", max.iter.harmony = 10)
sce2 <- sce2 %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = 0.02, verbose = FALSE)


sce2$group <- sce2$sub_atlas
sce2@meta.data[sce2$seurat_clusters=='1' & sce2$sub_atlas=='SCLC','group'] <- 'NSCLC-like SCLC'
sce2$group2 <- sce2$group
sce2@meta.data[sce2$group %in% c("LUAD","LUSC"),'group2'] <- 'NSCLC'

DimPlot(sce2,reduction = "umap",cols = c("#B24745FF","#79AF97FF","#7B4173FF"),
        label = F,group.by = "group2")
ggplot2::ggsave(filename = 'nsclc_sclc_umap.pdf',width = 8,height = 7)
saveRDS(sce2,file = 'nsclc_sclc.rds')
library(MuDataSeurat)
MuDataSeurat::WriteH5AD(sce2, 'nsclc_sclc.h5ad', assay="RNA")



## sclc marker visualization 
scmarker <- c("ASCL1","NEUROD1","POU2F3","YAP1")
VlnPlot(sce2,scmarker,group.by = 'group2',flip = F,stack = T,split.by = "group2")+scale_fill_jama()
ggsave("nsclc_sclc_sclcmaker.pdf")

## neu marker visualization
neuromarker <- c("DLL3","INSM1","HES6","CHGA")
VlnPlot(sce2,neuromarker,group.by = 'group2',flip = F,stack = T,split.by = "group2")+scale_fill_jama()
ggsave("nsclc_sclc_neuromaker.pdf")



# CytoTRACE for sclc
library(CytoTRACE)
library(ggprism)
library(ggpubr)

results <- CytoTRACE(as.matrix(sce2@assays$RNA@counts), ncores = 8)
dat <- data.frame(group=sce2$group2,cytotrace=results$CytoTRACE)
ggplot(dat,aes(x=group,y=cytotrace))+
  stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+
  geom_boxplot(aes(fill=group),
               outlier.colour="white",size=0.8)+
  theme(panel.background =element_blank(), 
        axis.line=element_line(),
        legend.position="none",plot.title = element_text(size=14))+
  ggtitle("boxplot")+
  stat_compare_means(comparisons = list(c("NSCLC", "NSCLC-like SCLC"),
                                        c("NSCLC-like SCLC","SCLC"),
                                        c("NSCLC","SCLC")),
                     size = 6,
                     method = "t.test",label = "p.value")+
  scale_fill_prism(palette = "candy_bright")
ggsave("nsclc_sclc_cyto_boxplot.pdf")

t.test(dat[dat$group=="NSCLC",2],dat[dat$group=="SCLC",2],var.equal=T)



## prepared for scarches mapping to spatial data
rawsce <- readRDS('epi_raw.rds')
sce <- readRDS('malignant_sub.rds')
subcds <- cds[,cds@clusters$UMAP$partitions == 1]
sclc_transition <- rawsce[,colnames(subcds)]
sclc_transition_tmp <- sce[,colnames(subcds)]
sclc_transition$stage <- sclc_transition_tmp$stage
sclc_transition <- sclc_transition[,sclc_transition$sub_atlas=='SCLC']
saveRDS(sclc_transition,'sclc_transition_symbol.rds')
library(MuDataSeurat)
MuDataSeurat::WriteH5AD(sclc_transition, 'sclc_transition_symbol.h5ad', assay="RNA")

library(clusterProfiler)
library(org.Hs.eg.db)
gene <- bitr(rownames(sclc_transition),
             fromType="SYMBOL",
             toType="ENSEMBL",
             OrgDb = "org.Hs.eg.db")
mat <- sclc_transition@assays$RNA$counts
mat2 <- mat[gene$SYMBOL,]
rownames(mat2) <- gene$ENSEMBL
sclc_transition2 <- CreateSeuratObject(counts = mat2, meta.data = sclc_transition@meta.data)
saveRDS(sclc_transition2,'sclc_transition_ensemble.rds')
library(MuDataSeurat)
MuDataSeurat::WriteH5AD(sclc_transition2, 'sclc_transition_ensemble.h5ad', assay="RNA")





