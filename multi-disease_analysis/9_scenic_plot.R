library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(SummarizedExperiment)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)



run_netplot <- function(celltype,regulon,ntop,node_size,lable_size){
  library(ggraph)
  grn <- fread(paste0(celltype,'_sub_latest_grn.tsv'),sep='\t',header=T,stringsAsFactors=F)
  inregulons <- gsub('[(+)]','',regulonsToPlot)
  c1 <- which(grn$TF %in% inregulons)
  grn2 <- grn[c1,]
  
  pdf(paste0(celltype,'_tf_netplot.pdf'),width = 10,height = 10)
  for (tf in unique(grn2$TF)) {
    tmp <- subset(grn2,TF==tf)
    if (dim(tmp)[1] > ntop) {
      tmp <- tmp[order(tmp$importance,decreasing=T),]
      tmp <- tmp[1:ntop,]
    }
    
    node2 <- data.frame(tmp$target)
    node2$node.size <- node_size
    node2$node.colour <- 'black'
    colnames(node2) <- c('node','node.size','node.colour')
    df1 <- data.frame(node=tf,node.size=2,node.colour='#A73030FF')
    node2 <- rbind(df1,node2)
    
    edge2 <- tmp
    colnames(edge2) <- c('from','to','edge.width')
    edge2$edge.colour <- "#5C615D"
    torange=c(0.1,1)
    edge2$edge.width <- scales::rescale(edge2$edge.width,to=torange)
    
    graph_data <- tidygraph::tbl_graph(nodes = node2, edges = edge2, directed = T)
    p1 <- ggraph(graph = graph_data, layout = "stress", circular = TRUE) + geom_edge_arc(aes(edge_colour = edge.colour, edge_width = edge.width)) +
      scale_edge_width_continuous(range = c(1,0.2)) +geom_node_point(aes(colour = node.colour, size = node.size))+ theme_void() +
      geom_node_label(aes(label = node,colour = node.colour),size = lable_size, repel = TRUE)
    p1 <- p1 + scale_color_manual(values=c('#A73030FF','black'))+scale_edge_color_manual(values=c("#5C615D"))
    print(p1)
  }
  dev.off()
}




load('mono_tf_scenic.rdata')
gse2 <- gse
rss <- calcRSS(AUC = getAUC(sub_regulonAUC), 
               cellAnnotation = cellTypes[colnames(sub_regulonAUC), 'mono_anno']) 
rss <- na.omit(rss) 
rssPlot <- plotRSS(rss)
regulonsToPlot <- c('KLF6(+)','KLF4(+)','XBP1(+)','FOS(+)','ETS1(+)')
regulonsToPlot %in% row.names(sub_regulonAUC)
gse2@meta.data = cbind(gse2@meta.data ,t(assay(sub_regulonAUC[regulonsToPlot,])))

# Visualization
VlnPlot(gse2, features = regulonsToPlot, group.by = 'mono_anno',
        pt.size = 0, cols = c("#895C48","#869CC3","#DEA164"))

FeaturePlot(gse2,features = regulonsToPlot, cols = c("grey", "#A73030FF"))


run_netplot(celltype = 'mono',
            regulon = regulonsToPlot,
            ntop = 50,
            node_size = 1.5,
            lable_size = 6)









load('bcell_tf_scenic.rdata')
gse2 <- gse
rss <- calcRSS(AUC = getAUC(sub_regulonAUC), 
               cellAnnotation = cellTypes[colnames(sub_regulonAUC), 'bcell_anno']) 
rss <- na.omit(rss) 
rssPlot <- plotRSS(rss)
regulonsToPlot <- c('POU2F2(+)','BCL11A(+)','BACH2(+)',
                    'TFDP1(+)','KLF6(+)',
                    'RELB(+)','FOXO1(+)')
regulonsToPlot %in% row.names(sub_regulonAUC)
gse2@meta.data = cbind(gse2@meta.data ,t(assay(sub_regulonAUC[regulonsToPlot,])))

# Visualization
DotPlot(gse2, features = regulonsToPlot, group.by = 'bcell_anno') +
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 90, hjust = 0.5,vjust=0.5))+
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))

FeaturePlot(gse2,features = regulonsToPlot, cols = c("grey", "#A73030FF"))


run_netplot(celltype = 'bcell',
            regulon = regulonsToPlot,
            ntop = 50,
            node_size = 1.5,
            lable_size = 6)




