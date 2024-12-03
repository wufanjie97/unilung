library(Seurat)
library(dplyr)
library(stringr)
library(tidyverse)
library(ggsci)


multi_barplot <- function(celltype){
  sce <- readRDS(paste0(celltype,'_sub_latest.rds'))
  
  ifelse(celltype=='mono',
         kk <- as.data.frame(table(sce$mono_anno,sce$donor_status)),
         kk <- as.data.frame(table(sce$bcell_anno,sce$donor_status)))

  ggplot(data = kk,aes(x=Var1,y=Freq,fill=Var2))+
    geom_bar(stat = "identity",
             position = "fill")+
    scale_fill_igv()+
    scale_y_continuous(expand = expansion(mult=c(0.01,0.1)),
                       labels = scales::percent_format())+
    theme(panel.background = element_blank(),
          axis.line = element_line(),
          legend.position = "right")+
    labs(x=NULL,y="Percent (%)")+
    guides(fill=guide_legend(title = NULL,ncol = 1,byrow = FALSE))
  ggplot2::ggsave(filename = file.path(dir,paste0(celltype,'_group_barplot.pdf')),width = 6,height = 10)
}

multi_barplot('mono')
multi_barplot('bcell')