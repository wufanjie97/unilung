library(Seurat)
library(dplyr)
library(EnhancedVolcano)
library(stringr)
library(tidyverse)

run_edgeR <- function(data, method = 'TMM') {
  library(Libra)
  library(MAST)
  library(edgeR)
  library(readxl)
  library(ggvenn)
  mat <- to_pseudobulk(data, cell_type_col = "second_anno", label_col = "deg_group", replicate_col = 'donor_ID')
  group <- data@meta.data %>% distinct(donor_ID, deg_group, .keep_all = T) %>% 
    dplyr::select(donor_ID, deg_group) %>% remove_rownames()
  rownames(group) <- paste0(group$donor_ID,':',group$deg_group)
  
  count <- mat[[1]][rownames(group)]
  group_list = group[match(colnames(count),rownames(group)),'deg_group']
  
  dge <- DGEList(counts = count, group = group_list) %>%
    calcNormFactors(method = method)
  
  design <- model.matrix(~0+group_list)
  rownames(design)<-colnames(dge)
  colnames(design)<-levels(group_list)
  
  dge <- estimateDisp(dge,design)
  fit <- glmFit(dge, design)
  fit2 <- glmLRT(fit, contrast=c(-1,1)) #-1: control; 1: experiment。
  
  res <- topTags(fit2, n=nrow(count)) %>% as.data.frame()
  return(res)
}


run_go <- function(deg,celltype,width=NULL,height=NULL){
  library(patchwork)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  deg_up <- subset(deg, logFC>1 & FDR<0.05)
  go_up <- enrichGO(gene = rownames(deg_up),
                    OrgDb = 'org.Hs.eg.db',
                    keyType = 'SYMBOL',
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.01,
                    qvalueCutoff = 0.01)
  
  go_top10 <- data.frame(go_up)[1:8,]
  go_top10$type <- "up"
  go_top10$qvalue <- -log10(go_top10$qvalue)
  ggplot(go_top10,aes(qvalue,Description))+
    geom_bar(aes(y=reorder(Description,qvalue),x=qvalue,fill=go_top10$type)
             ,stat='identity',width=0.5)+
    scale_fill_manual(values = c(down="#7DC5A0",up="#D58890"))+
    theme_bw()+
    theme(panel.grid = element_blank())
  ggplot2::ggsave(filename = paste0(celltype,"_go.pdf"),width = width,height = height)
}

# DEG analysis between lym-monocyte and others
sce2 <- readRDS('mono_sub_latest.rds')

Idents(sce2) <- sce2$mono_anno
sce2$deg_group <- ifelse(sce2$mono_anno == 'Lym-monocyte', 'Lym-monocyte', "Others")
sce2$deg_group <- factor(sce2$deg_group, levels = c('Others', 'Lym-monocyte'))

mono_deg_bulk <- run_edgeR(data = sce2, method = 'TMM')

Idents(sce2) <- sce2$deg_group
mono_deg_sc <- FindMarkers(sce2,ident.1 = 'Lym-monocyte',
                           ident.2 = 'Others',logfc.threshold = 0.1,test.use = "wilcox")


deg_list <- mono_deg_bulk[intersect(rownames(mono_deg_bulk),rownames(mono_deg_sc)),]
EnhancedVolcano(deg_list,
                lab = rownames(deg_list),
                selectLab = c('TMSB4XP6','ZEB1','LIME1','TRAF1'),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.00001,
                pointSize = 2,
                labSize = 4,
                xlim = c(-2, 2),
                ylim = c(0,200),
                col = c('#b5b5b5','#b5b5b5','#4D4398','#F18D00'),
                gridlines.major = FALSE,
                gridlines.minor = FALSE)


## GO enrichment
run_go(mono_deg_bulk,'mono',width = 10,height=8)




# DEG analysis between T-like B cell and others
sce2 <- readRDS('bcell_sub_latest.rds')

Idents(sce2) <- sce2$bcell_anno
sce2$deg_group <- ifelse(sce2$bcell_anno == 'T-like B cell', 'T-like B cell', "Others")
sce2$deg_group <- factor(sce2$deg_group, levels = c('Others', 'T-like B cell'))

bcell_deg_bulk <- run_edgeR(data = sce2, method = 'TMM')

Idents(sce2) <- sce2$deg_group
bcell_deg_sc <- FindMarkers(sce2,ident.1 = 'T-like B cell',
                           ident.2 = 'Others',logfc.threshold = 0.1,test.use = "wilcox")




deg_list <- bcell_deg_bulk[intersect(rownames(bcell_deg_bulk),rownames(bcell_deg_sc)),]
EnhancedVolcano(deg_list,
                lab = rownames(deg_list),
                selectLab = c('TRBC1','GZMA','RGS1','HLA-DQB1','CD83'),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.00001,
                pointSize = 2,
                labSize = 4,
                xlim = c(-2, 2),
                ylim = c(0,400),
                col = c('#b5b5b5','#b5b5b5','#4D4398','#F18D00'),
                gridlines.major = FALSE,
                gridlines.minor = FALSE)


## GO enrichment
run_go(bcell_deg_bulk,'bcell',width = 10,height=8)



