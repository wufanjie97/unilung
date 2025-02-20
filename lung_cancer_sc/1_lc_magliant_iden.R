library(reticulate)
library(Seurat)
library(tidyverse)
library(harmony)
library(cowplot)
dir <- "lung_cancer_sc"
setwd(dir)

use_condaenv("envs/scvi")
scanpy <- import("scanpy")
pd <- import("pandas")
np <- import("numpy")


in_ad <- scanpy$read_h5ad(file.path(dir,"epi_raw.h5ad"))
count <- np$transpose(in_ad$X)
colnames(count) <- rownames(in_ad$obs)
rownames(count) <- rownames(in_ad$var)
seu <- CreateSeuratObject(counts = count, meta.data = in_ad$obs)
saveRDS(seu,file.path(dir,"epi_raw.rds"))



sce <- NormalizeData(seu,normalization.method = "LogNormalize",scale.factor = 10000)
sce <- FindVariableFeatures(sce,reductionselection.method = "vst",nfeatures = 3000)
sce <- ScaleData(sce)
sce <- RunPCA(object = sce,npcs = 30,pc.genes=VariableFeatures(sce),verbose = F)
sce <- RunHarmony(sce, group.by.vars = "donor_ID", max.iter.harmony = 10)
sce <- sce %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = 0.1)



gse_marker <- RenameIdents(object = sce,
                           "0" = "Type II alveolar cell", "1" = "Basal cell",
                           "2" = "Malignant cell", "3" = "Malignant cell",
                           "4" = "Malignant cell", "5" = "Malignant cell",
                           "6" = "Ciliated cell", "7" = "Type I alveolar cell"
)
sce$annotation<-Idents(gse_marker)


select_genes <- c(
  "KRT13","KRT17","S100A2",#Basal cell
  "CAV1","AGER",#Alveolar type 1
  "SFTPC","SFTPA1","SFTPA2",#Alveolar type 2
  "SCGB1A1","SCGB3A1",#club cell
  "TPPP3","FOXJ1","PIFO",#ciliated cell
  "EPCAM"# cancer
)


DotPlot(sce,features = select_genes,group.by = "annotation")+
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))


library(ggsci)
library(cowplot)
plot_grid(ncol=2,
          DimPlot(sce2,reduction = "umap",cols = pal_nejm("default")(6),
                  label = F,group.by = "annotation")+ggtitle("Cell type"),
          DimPlot(sce2,reduction = "umap",cols = c("#B24745FF","#79AF97FF","#7B4173FF"),
                  label = F,group.by = "sub_atlas")+ggtitle("Cancer type"))


saveRDS(sce,file.path(dir,"epi.rds"))
library(MuDataSeurat)
MuDataSeurat::WriteH5AD(sce, 'epi.h5ad', assay="RNA")





## infercnv
library(infercnv)
library(AnnoProbe)
library(reticulate)
dir <- "lung_cancer_sc"

use_condaenv("/envs/scvi")
scanpy <- import("scanpy")
np <- import("numpy")

sce <- readRDS(file.path(dir,'epi.rds'))
ref_epi <- scanpy$read_h5ad(file.path(dir,"ref_epi10k.h5ad"))
epi_dat <- np$transpose(ref_epi$X)
colnames(epi_dat) <- rownames(ref_epi$obs)
rownames(epi_dat) <- rownames(ref_epi$var)
save(sce,epi_dat,file = file.path(dir,"infercnv/epidata.RData"))



run_cnv <- function(times,cluster){
  library(Seurat)
  library(infercnv)
  library(AnnoProbe)
  library(reticulate)
  options(future.globals.maxSize= 5*1024^3)
  options(scipen = 100)
  
  load(file.path(dir,"infercnv/epidata.RData"))
  gse <- sce[,sce@meta.data$seurat_clusters %in% cluster]
  dat <- as.data.frame(gse@assays$RNA@counts)
  colnames(dat) <- colnames(gse)
  rownames(dat) <- rownames(gse)
  
  dat2 <- dat[intersect(rownames(dat),rownames(epi_dat)),]
  epi_dat2 <- epi_dat[intersect(rownames(dat),rownames(epi_dat)),]
  
  all_dat <- cbind(dat2,epi_dat2)
  groupinfo <- data.frame(v1=colnames(all_dat),
                          v2=c(gse$seurat_clusters,rep('ref-epithelial',10000)))
  
  geneInfor <- annoGene(rownames(all_dat),"SYMBOL",'human')
  geneInfor <- geneInfor[with(geneInfor,order(chr,start)),c(1,4:6)]
  geneInfor <- geneInfor[!duplicated(geneInfor[,1]),]
  
  data <- all_dat[rownames(all_dat) %in% geneInfor[,1],]
  data <- data[match(geneInfor[,1], rownames(data)),] 
  
  print("=============== Start save files ===============")
  expfile <- file.path(paste0(dir,'/infercnv'),paste0('expfile',times,'.txt'))
  write.table(data,file = expfile,sep = '\t',quote = F)
  groupfile <- file.path(paste0(dir,'/infercnv'),paste0('groupfile',times,'.txt'))
  write.table(groupinfo,file = groupfile,sep = '\t',quote = F,col.names = F,row.names = F)
  genefile <- file.path(paste0(dir,'/infercnv'),paste0('genefile',times,'.txt'))
  write.table(geneInfor,file = genefile,sep = '\t',quote = F,col.names = F,row.names = F)
  print("=============== Done save files ===============")

  
  options(scipen = 100)
  options("Seurat.object.assay.version" = "v3")
  
  expfile <- file.path(paste0(dir,'/infercnv'),paste0('expfile',times,'.txt'))
  groupfile <- file.path(paste0(dir,'/infercnv'),paste0('groupfile',times,'.txt'))
  genefile <- file.path(paste0(dir,'/infercnv'),paste0('genefile',times,'.txt'))

  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=expfile,
                                      annotations_file=groupfile,
                                      delim="\t",
                                      gene_order_file= genefile,
                                      ref_group_names='ref-epithelial') 
  print("=============== Done create cnvdata ===============")
  future::plan("multicore",workers=16)
  infercnv_obj2 = infercnv::run(infercnv_obj,
                                cutoff=0.1, 
                                out_dir=paste0(dir,"/infercnv/cnv",times),
                                cluster_by_groups=T,
                                write_expr_matrix=T,
                                write_phylo = T, 
                                denoise = T,
                                HMM = F,
                                output_format = "pdf")
  save(infercnv_obj,infercnv_obj2,file = paste0(dir,"/infercnv/infercnv_obj",times,".RData"))
}

## Due to the scale of data, We split it into three subsets according to the clustering results and run them separately
run_cnv(times = 1, cluster = c(0,16:34))
run_cnv(times = 2, cluster = c(1:3,13:15))
run_cnv(times = 3, cluster = c(4:12))






library(infercnv)
library(Seurat)
library(ggplot2)
library(data.table)
library(tidyverse)
dir <- "lung_cancer_sc"
setwd("lung_cancer_sc/infercnv")


cnvScore <- function(data){
  data <- data %>% as.matrix() %>%
    t() %>% 
    scale() %>% 
    rescale(to=c(-1, 1)) %>% 
    t()
  
  cnv_score <- as.data.frame(colSums(data * data))
  return(cnv_score)
}

# dat <- read.table("infercnv.observations.txt", header=T)
dat <- readRDS("cnv_all/run.final.infercnv_obj")
cnv_score <- cnvScore(dat@expr.data)

colnames(cnv_score) <- "score"
cnv_score <- rownames_to_column(cnv_score, var='cellid')

group <- read.table("groupfile.txt")
colnames(group) <- c("cellid","cluster")

cnv_score <- left_join(cnv_score,group,by="cellid")


ggplot(cnv_score,aes(x=cluster, y=score,fill=cluster))+geom_violin(aes(fill=cluster),color="NA")+
  scale_fill_igv()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",linewidth = .4),
        axis.ticks.x=element_line(color="black",linewidth = .4,lineend = 1),
        axis.ticks.y=element_line(color="black",linewidth = .4,lineend = 1),
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 9))+
  ylab("cnv score")+xlab("celltype")
ggsave("cnv_score.pdf",width = 15,height = 6,units = "cm")




library(CytoTRACE)
dir <- "lung_cancer_sc"
setwd(dir)


sce <- readRDS(file.path(dir,'epi.rds'))
sce <- FindVariableFeatures(sce,reductionselection.method = "vst",nfeatures = 2000)
sce2 <- SketchData(object = sce,ncells = 100000,method = "LeverageScore",sketched.assay = "bpcell")
sce2 <- sce2@assays$bpcell@counts
sce_sub <- subset(sce,cells=colnames(sce2))

mat <- as.matrix(sce_sub@assays$RNA@counts)
results <- CytoTRACE(mat, ncores = 8)

emb <- sce_sub@reductions[["umap"]]@cell.embeddings

celltype <- as.character(sce_sub$annotation)
names(celltype) <- rownames(sce_sub@meta.data)
plotCytoTRACE(results, phenotype = celltype,emb = emb,outputDir = paste0(dir,"/harm_cytotrace/anno_"))

atlas <- as.character(sce_sub$sub_atlas)
names(atlas) <- rownames(sce_sub@meta.data)
plotCytoTRACE(results, phenotype = atlas, emb = emb, outputDir = paste0(dir,"/harm_cytotrace/atlas_"))


