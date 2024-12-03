# devtools::install_github("liuqivandy/SpaGene")
library(Seurat)
library(tidyverse)
library(SpaGene)
library(cowplot)
options(Seurat.object.assay.version = "v3")
dir <- 'lung_cancer_sp'
options(future.globals.maxSize = 10000 * 1024^2)
setwd(dir)



rds_from_h5ad <- function(dataset){
  library(reticulate)
  use_condaenv("envs/cell2loc")
  scanpy <- import("scanpy")
  pd <- import("pandas")
  np <- import("numpy")
  cancer <- ifelse(grepl('luad', dataset), 'LUAD', 'LUSC')
  for (sample in list.files(paste0(dataset, '/E-MTAB-13530'))) {
    data <- Load10X_Spatial(data.dir = paste0(dataset, '/E-MTAB-13530/', sample, '/'),
                            filename = "filtered_feature_bc_matrix.h5",
                            assay = "Spatial",
                            slice = basename(sample))
    ## prepare python enviroment
    use_condaenv("envs/cell2loc")
    scanpy <- import("scanpy")
    pd <- import("pandas")
    np <- import("numpy")
    ## load h5ad
    in_ad <- scanpy$read_h5ad(paste0('spagene/h5adobj/',dataset,'_',sample,'.h5ad'))
    meta <- in_ad$obs
    ## the barcode in h5ad is raw barcode+'-1' (eg.AAACATTTCCCGGATT-1-1)
    rownames(meta) <- gsub("(-[^-]*)$", "", gsub("(-[^-]*)", "\\1", rownames(meta)))
    data2 <- data[,rownames(meta)]
    data2@meta.data <- cbind(data2@meta.data,meta)
    saveRDS(data2,file = paste0('spagene/rdsobj/',dataset,'_',sample,'.rds'))
  }
}

rds_from_h5ad('luad')
rds_from_h5ad('lusc')






## NSCLC-like SCLC region share gene module in luad and lusc
lusc1 <- readRDS('lusc_P17_T2.rds')
lusc2 <- readRDS('lusc_P19_T1.rds')
luad1 <- readRDS('luad_P10_T1.rds')
luad2 <- readRDS('luad_P24_T2.rds')

lusc1 <- readRDS('spagene/rdsobj/lusc_P17_T2.rds')
lusc2 <- readRDS('spagene/rdsobj/lusc_P19_T1.rds')
luad1 <- readRDS('spagene/rdsobj/luad_P10_T1.rds')
luad2 <- readRDS('spagene/rdsobj/luad_P24_T2.rds')

count1 <- GetAssayData(lusc1,slot="counts")
count2 <- GetAssayData(lusc2,slot="counts")
count3 <- GetAssayData(luad1,slot="counts")
count4 <- GetAssayData(luad2,slot="counts")


location1 <- GetTissueCoordinates(lusc1)
location2 <- GetTissueCoordinates(lusc2)
location3 <- GetTissueCoordinates(luad1)
location4 <- GetTissueCoordinates(luad2)

spa1 <- SpaGene(count1,location1)
spa2 <- SpaGene(count2,location2)
spa3 <- SpaGene(count3,location3)
spa4 <- SpaGene(count4,location4)

pattern <- FindPattern_Multi(list(spa1,spa2,spa3,spa4),nPattern=10)
locationlist <- list(location1[,2:1],location2[,2:1],location3[,2:1],location4[,2:1])




## choose pattern
pdf("spagene/share_pattern.pdf", width = 8, height = 8)
PlotPattern_Multi(pattern,locationlist,pt.size=1,patternid=4,max.cutoff =0.95)
PlotPattern_Multi(pattern,locationlist,pt.size=1,patternid=6,max.cutoff =0.95)
PlotPattern_Multi(pattern,locationlist,pt.size=1,patternid=7,max.cutoff =0.95)
dev.off()


top5 <- apply(pattern$genepattern,2,function(x){names(x)[order(x,decreasing=T)][1:5]})
library(pheatmap)
pheatmap(pattern$genepattern[rownames(pattern$genepattern)%in%top5,],
         fontsize_row = 6,
         cellwidth = 8, 
         cellheight = 10, 
         fontsize = 8, 
         cluster_col = F,
         filename = "spagene/share_patterns_heatmap.pdf")



save(spa1,spa2,spa3,spa4,pattern,file = 'spagene/pattern.rdata')




## identify colocalized ligand-receptor pairs
load("spagene/LRpair_human.rds")
lr1 <- SpaGene_LR(count1,location1,LRpair=LRpair)
write.csv(lr1,'lusc_P17_T2_LRI.csv',row.names = T)
lr2 <- SpaGene_LR(count2,location2,LRpair=LRpair)
write.csv(lr2,'lusc_P19_T1_LRI.csv',row.names = T)
lr3 <- SpaGene_LR(count3,location3,LRpair=LRpair)
write.csv(lr3,'luad_P10_T1_LRI.csv',row.names = T)
lr4 <- SpaGene_LR(count4,location4,LRpair=LRpair)
write.csv(lr4,'luad_P24_T2_LRI.csv',row.names = T)

save(count1,count2,count3,count4,
     location1,location2,location3,location4,
     lr1,lr2,lr3,lr4,LRpair,file = 'spagene/LRI.rdata')


get_top_rows <- function(df,ntop) {
  top_rows <- df[order(df$adj), ][1:ntop, ]
  return(rownames(top_rows))
}

ntop = 30
top_rows_lr1 <- get_top_rows(lr1,ntop)
top_rows_lr2 <- get_top_rows(lr2,ntop)
top_rows_lr3 <- get_top_rows(lr3,ntop)
top_rows_lr4 <- get_top_rows(lr4,ntop)

# shared LRI
share_LRI <- Reduce(intersect, list(top_rows_lr1, top_rows_lr2, top_rows_lr3, top_rows_lr4))


# choose LRI "COL1A1_ITGB1" "COL1A2_ITGB1" "FN1_ITGB1" "LUM_ITGB1"
library(cowplot)
for (i in c("COL1A1","COL1A2","FN1","LUM")) {
  plot_grid(ncol=2,
            plotLR(count1,location1,LRpair=c(i,"ITGB1"),alpha.min=0.2,pt.size = 1.5),
            plotLR(count2,location2,LRpair=c(i,"ITGB1"),alpha.min=0.2,pt.size = 1.5),
            plotLR(count3,location3,LRpair=c(i,"ITGB1"),alpha.min=0.2,pt.size = 1.5),
            plotLR(count4,location4,LRpair=c(i,"ITGB1"),alpha.min=0.2,pt.size = 1.5))
  ggplot2::ggsave(filename = paste0('spagene/share_LRI_',i,'.pdf'),width = 15,height = 7)
}


plotLR(count1,location1,LRpair=c("COL1A1","ITGB1"),alpha.min=0.2,pt.size = 1.5)

