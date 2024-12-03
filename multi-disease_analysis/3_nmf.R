library(Seurat)
library(tidyverse)
library(cowplot)
library(NMF)

dir <- "multi-disease_analysis"

nmfpatterns <- function(celltype){
  library(NMF)
  library(ggplot2)
  library(Seurat)
  library(dplyr)
  library(viridis)
  ranks <- 2:8
  sce <- readRDS(file.path(dir,paste0(celltype,"_sub_raw.rds")))
  
  sce <- FindVariableFeatures(sce,selection.method = "vst",nfeatures = 2000)
  sce <- ScaleData(sce)
  sce <- RunPCA(object = sce,npcs = 30,pc.genes=VariableFeatures(sce),verbose = F)
  sce_hvg <- sce[VariableFeatures(sce),]
  
  nmf_counts <- GetAssayData(object = sce_hvg, slot = "counts") %>% as.matrix()
  
  # remove the columns that are all zeros or have missing values
  i0 <- which(colSums(nmf_counts) == 0)
  i_na <- which(colSums(is.na(nmf_counts)) > 0)
  # remove the row which contains at least one null or NA-filled
  if (length(i0)>0 | length(i_na)>0) {
    nmf_counts <- nmf_counts[which(rowSums(nmf_counts) > 0),-c(i0, i_na)] %>% log1p()
  } else {
    nmf_counts <- nmf_counts[which(rowSums(nmf_counts) > 0),] %>% log1p()
  }
  
  cophenetic_coeffs <- numeric(length(ranks))
  
  for (i in seq_along(ranks)) {
    rank <- ranks[i]
    nmf_results <- nmf(nmf_counts, rank=rank, nrun = 20, .options = 'vP4', .pbackend = 4, seed = 123456)
    cophenetic_coeffs[i] <- cophcor(nmf_results)
  }
  cophenetic_df <- data.frame(FactorizationRank = ranks, CopheneticCoefficient = cophenetic_coeffs)
  save(nmf_results,cophenetic_df,cophenetic_coeffs,file = file.path(dir,paste0(celltype,"_nmf_chooserank.RData")))
  
  ggplot(cophenetic_df, aes(x = FactorizationRank, y = CopheneticCoefficient)) +
    geom_point(color = "purple", size = 3) +
    geom_line(color = "purple") +
    labs(title = "Cophenetic correlation", x = "Factorization rank", y = "Coefficient") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 16))
  ggplot2::ggsave(filename = file.path(dir,paste0(celltype,"_NMFpattern.pdf")),width = 10,height = 7)
  return(cophenetic_df)
}

mono_nmf_df <- nmfpatterns("mono")
bcell_nmf_df <- nmfpatterns("bcell")



## NMF for monocyte
celltype='mono'
sce <- readRDS(file.path(dir,paste0(celltype,"_sub_raw.rds")))
nmf_counts <- GetAssayData(object = sce, slot = "counts") %>% as.matrix()
i0 <- which(colSums(nmf_counts) == 0)
i_na <- which(colSums(is.na(nmf_counts)) > 0)
if (length(i0)>0 | length(i_na)>0) {
  nmf_counts <- nmf_counts[which(rowSums(nmf_counts) > 0),-c(i0, i_na)] %>% log1p()
} else {
  nmf_counts <- nmf_counts[which(rowSums(nmf_counts) > 0),] %>% log1p()
}
nmf_results <- nmf(nmf_counts, rank=3, nrun = 20, .options = 'vP4', .pbackend = 4, seed = 123456)
save(nmf_counts,nmf_results,file = paste0(celltype,'_nmf.rdata'))


## bcell
celltype='bcell'
sce <- readRDS(file.path(dir,paste0(celltype,"_sub_raw.rds")))
nmf_counts <- GetAssayData(object = sce, slot = "counts") %>% as.matrix()
i0 <- which(colSums(nmf_counts) == 0)
i_na <- which(colSums(is.na(nmf_counts)) > 0)
if (length(i0)>0 | length(i_na)>0) {
  nmf_counts <- nmf_counts[which(rowSums(nmf_counts) > 0),-c(i0, i_na)] %>% log1p()
} else {
  nmf_counts <- nmf_counts[which(rowSums(nmf_counts) > 0),] %>% log1p()
}
nmf_results <- nmf(nmf_counts, rank=4, nrun = 20, .options = 'vP4', .pbackend = 4, seed = 123456)
save(nmf_counts,nmf_results,file = paste0(celltype,'_nmf.rdata'))
