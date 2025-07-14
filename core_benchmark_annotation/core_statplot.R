library(ggplot2)
library(dplyr)
library(tidyverse)
library(stringr)
library(ggforce)


meta <- read.csv("core_meta.csv",row.names = 1)



## disease bar plot
status <- as.data.frame(table(meta$donor_status))
status$per <- status$Freq / sum(status$Freq) * 100
colnames(status) <- c("disease","num","per")
status$disease <- as.character(status$disease)


ggplot(status,aes(x,per))+
  geom_bar(stat="identity",width=0.4,aes(y=per,x=reorder(disease, -per)))+
  theme_bw()+
  theme(axis.text.x = element_text(colour = 'black',size = 10, angle = 60, hjust = 1, vjust = 1),
        axis.text.y = element_text(colour = 'black',size = 10),
        axis.title.x=element_text(size=12, colour = 'black'),
        axis.title.y=element_text(size=12, colour = 'black'),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))+
  xlab(label = 'Status')+ylab(label = 'Percent of cells')+
  guides(fill=FALSE)
ggplot2::ggsave("core_stat_disease.pdf",width = 3, height = 3)


## seq_tech bar plot
seq <- as.data.frame(table(meta$seq_tech))
seq$per <- seq$Freq / sum(seq$Freq) * 100
colnames(seq) <- c("seq","num","per")
seq$seq <- as.character(seq$seq)


ggplot(seq)+
  geom_bar(stat="identity",width=.5,aes(y=per,x=seq))+
  theme_bw()+
  theme(axis.text.x = element_text(colour = 'black',size = 10, angle = 60, hjust = 1, vjust = 1),
        axis.text.y = element_text(colour = 'black',size = 10),
        axis.title.x=element_text(size=12, colour = 'black'),
        axis.title.y=element_text(size=12, colour = 'black'),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))+
  xlab(label = 'Platform')+ylab(label = 'Percent of cells')+
  guides(fill=FALSE)
ggplot2::ggsave("core_stat_seq.pdf",width = 3, height = 4)


## age bar plot
age <- data.frame(age=c("0-10 year","11-20 year","21-30 year","31-40 year",
                        "41-50 year","51-60 year","61-70 year","71-80 year","81+ year","Unclassified"),
                  num=c(27684,3077,24682,31876,34973,101953,221677,117836,13608,277569))
age$per <- age$num / sum(age$num) * 100
age$x <- factor(age$age,levels=c("0-10 year","11-20 year","21-30 year","31-40 year",
                                 "41-50 year","51-60 year","61-70 year","71-80 year","81+ year","Unclassified"))


ggplot(age)+
  geom_bar(stat="identity",width=.5,aes(y=per,x=x))+
  theme_bw()+
  theme(axis.text.x = element_text(colour = 'black',size = 10, angle = 60, hjust = 1, vjust = 1),
        axis.text.y = element_text(colour = 'black',size = 10),
        axis.title.x=element_text(size=12, colour = 'black'),
        axis.title.y=element_text(size=12, colour = 'black'),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))+
  xlab(label = 'Age')+ylab(label = 'Percent of cells')+
  guides(fill=FALSE)
ggplot2::ggsave("core_stat_age.pdf",width = 3, height = 3)


## ethnicity bar plot
race <- as.data.frame(table(meta$ethnicity))
race$per <- race$Freq / sum(race$Freq) * 100
colnames(race) <- c("eth","num","per")
race$eth <- as.character(race$eth)


ggplot(race)+
  geom_bar(stat="identity",width=.5,aes(y=per,x=eth))+
  theme_bw()+
  theme(axis.text.x = element_text(colour = 'black',size = 10, angle = 60, hjust = 1, vjust = 1),
        axis.text.y = element_text(colour = 'black',size = 10),
        axis.title.x=element_text(size=12, colour = 'black'),
        axis.title.y=element_text(size=12, colour = 'black'),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))+
  xlab(label = 'Ethnicity')+ylab(label = 'Percent of cells')+
  guides(fill=FALSE)
ggplot2::ggsave("core_stat_ethnic.pdf",width = 2, height = 3)



## gender bar plot
gender <- as.data.frame(table(meta$donor_gender))
gender$per <- gender$Freq / sum(gender$Freq) * 100
colnames(gender) <- c("gen","num","per")
gender$gen <- as.character(gender$gen)
gender$gen <- c("Female","Male","Unclassified")


ggplot(gender)+
  geom_bar(stat="identity",width=.5,aes(y=per,x=gen))+
  theme_bw()+
  theme(axis.text.x = element_text(colour = 'black',size = 10, angle = 60, hjust = 1, vjust = 1),
        axis.text.y = element_text(colour = 'black',size = 10),
        axis.title.x=element_text(size=12, colour = 'black'),
        axis.title.y=element_text(size=12, colour = 'black'),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))+
  xlab(label = 'Sex')+ylab(label = 'Percent of cells')+
  guides(fill=FALSE)
ggplot2::ggsave("core_stat_sex.pdf",width = 1.5, height = 3)





## all celltype bar plot
celltype <- as.data.frame(table(meta$unilung_final_celltype))
celltype$per <- celltype$Freq / sum(celltype$Freq) * 100
colnames(celltype) <- c("celltype","num","per")
celltype$celltype <- as.character(celltype$celltype)


ggplot(celltype)+
  geom_bar(stat="identity",width=.5,aes(y=num,x=reorder(celltype, -num)))+
  theme_bw()+
  theme(axis.text.x = element_text(colour = 'black',size = 10, angle = 60, hjust = 1, vjust = 1),
        axis.text.y = element_text(colour = 'black',size = 10),
        axis.title.x=element_text(size=10, colour = 'black'),
        axis.title.y=element_text(size=10, colour = 'black'),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))+
  xlab(label = 'Celltype')+ylab(label = 'Number of cells')+
  guides(fill=FALSE)
ggplot2::ggsave("core_stat_celltype.pdf",width = 7, height = 4)






### sanky plot for cell annotation level1
library(ggplot2)
library(ggalluvial)

uhaf_cols <- c("unilung_ann_level1", "unilung_ann_level2", "unilung_ann_level3", "unilung_ann_level4",
               "unilung_ann_level2_final", "unilung_ann_level3_final", "unilung_ann_level4_final")
map_unilung_cols <- c("map_ann_level1", "unilung_ann_level1")

if (!all(uhaf_cols %in% colnames(meta))) stop("Some uhaf columns missing in meta.")
if (!all(map_unilung_cols %in% colnames(meta))) stop("Some map_unilung columns missing in meta.")

uhaf <- meta[, uhaf_cols]
uhaf$unilung_ann_level2 <- gsub("Unclassified", "", uhaf$unilung_ann_level2)
uhaf$unilung_ann_level3 <- gsub("Unclassified", "", uhaf$unilung_ann_level3)
uhaf$unilung_ann_level4 <- gsub("Unclassified", "", uhaf$unilung_ann_level4)

map_unilung <- meta[, map_unilung_cols]

meta$cell_id <- 1:nrow(meta) 
df_uhaf <- to_lodes_form(uhaf, axes = 1:4, id = "cell_id")
df_map_unilung <- to_lodes_form(map_unilung, axes = 1:2, id = "cell_id")
df <- df_map_unilung

col <- rep(c('#DF605E', '#9EBB65', '#63B8B8', '#5579B5', '#F6B56A', '#794B21', 
             '#DCBD9B', '#E8CC32', '#BE752D', '#4a8594'), 5)


pdf("map_unilung_sankey_lv1.pdf",width = 8, height = 6)

ggplot(df, aes(x = x, stratum = stratum, alluvium = cell_id, fill = stratum)) +
  geom_flow(width = 0.3, curve_type = "sine", alpha = 0.5, color = 'white', size = 0.1) +
  geom_stratum(width = 0.05) +
  geom_text(aes(label = stratum), stat = 'stratum', size = 2, color = 'black') + 
  scale_fill_manual(values = col) +
  theme_void() +
  theme(legend.position = 'none') +
  labs(title = "Sankey Diagram: map_ann_level1 to unilung_ann_level1")

dev.off()




### sanky plot for cell annotation level2
uhaf_cols <- c("unilung_ann_level1", "unilung_ann_level2", "unilung_ann_level3", "unilung_ann_level4",
               "unilung_ann_level2_final", "unilung_ann_level3_final", "unilung_ann_level4_final")
map_unilung_cols <- c("map_ann_level2", "unilung_ann_level2_final")


if (!all(uhaf_cols %in% colnames(meta))) stop("Some uhaf columns missing in meta.")
if (!all(map_unilung_cols %in% colnames(meta))) stop("Some map_unilung columns missing in meta.")


meta_filter <- meta[meta$map_ann_level2 != "Unclassified", ]
uhaf <- meta_filter[, uhaf_cols]
uhaf$unilung_ann_level2 <- gsub("Unclassified", "", uhaf$unilung_ann_level2)
uhaf$unilung_ann_level3 <- gsub("Unclassified", "", uhaf$unilung_ann_level3)
uhaf$unilung_ann_level4 <- gsub("Unclassified", "", uhaf$unilung_ann_level4)

map_unilung <- meta_filter[, map_unilung_cols]


meta_filter$cell_id <- 1:nrow(meta_filter) 
df_uhaf <- to_lodes_form(uhaf, axes = 1:4, id = "cell_id")
df_map_unilung <- to_lodes_form(map_unilung, axes = 1:2, id = "cell_id")


df <- df_map_unilung
col <- rep(c('#DF605E', '#9EBB65', '#63B8B8', '#5579B5', '#F6B56A', '#794B21', 
             '#DCBD9B', '#E8CC32', '#BE752D', '#4a8594'), 5)


pdf("map_unilung_sankey_lv2.pdf",width = 6, height = 10)

ggplot(df, aes(x = x, stratum = stratum, alluvium = cell_id, fill = stratum)) +
  geom_flow(width = 0.3, curve_type = "sine", alpha = 0.5, color = 'white', size = 0.1) +
  geom_stratum(width = 0.05,color = NA) +
  geom_text(aes(label = stratum), stat = 'stratum', size = 2, color = 'black') + 
  scale_fill_manual(values = col) +
  theme_void() +
  theme(legend.position = 'none') +
  labs(title = "Sankey Diagram: map_ann_level1 to unilung_ann_level1")

dev.off()



