library(ggplot2)
library(reshape2)
library(ggsci)


## bar plot
perlv1 <- read.csv('anno_stat_lv1.csv')
perlv2 <- read.csv('anno_stat_lv2.csv')
perlv3 <- read.csv('anno_stat_lv3.csv')
perlv4 <- read.csv('anno_stat_lv4.csv')


data = melt(perlv1)
colnames(data) = c('celltype','labelled','percent')

ggplot(data, aes(x = celltype, weight = percent, fill = labelled))+
  geom_bar(position = "stack",width=.4)+
  scale_fill_manual(values=c(
    match_labelled = "#BC3C29FF", cover_labelled = "#E18727FF", unmatch_labelled = "#0072B5FF"
    ))+
  # scale_fill_nejm()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        plot.title=element_text(color="black", size=10, vjust=0.5),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",linewidth = .4),
        axis.ticks.x=element_line(color="black",linewidth = .4,lineend = 1),
        axis.ticks.y=element_line(color="black",linewidth = .4,lineend = 1),
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 9))+
  labs(x="unilung annotation",y="Percentage of cells",title="Comparison between original and unilung label (level1)")

data = perlv2[-27,]
data = melt(data)
colnames(data) = c('celltype','labelled','percent')


ggplot(data, aes(x = celltype, weight = percent, fill = labelled))+
  geom_bar(position = "stack",width=.4)+
  scale_fill_manual(values=c(
    match_labelled = "#BC3C29FF", cover_labelled = "#E18727FF", unmatch_labelled = "#0072B5FF"
  ))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        plot.title=element_text(color="black", size=10, vjust=0.5),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",linewidth = .4),
        axis.ticks.x=element_line(color="black",linewidth = .4,lineend = 1),
        axis.ticks.y=element_line(color="black",linewidth = .4,lineend = 1),
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 9))+
  labs(x="unilung annotation",y="Percentage of cells",title="Comparison between original and unilung label (level2)")


data = perlv3[-26,]
data = melt(data)
colnames(data) = c('celltype','labelled','percent')

ggplot(data, aes(x = celltype, weight = percent, fill = labelled))+
  geom_bar(position = "stack",width=.4)+
  scale_fill_manual(values=c(
    match_labelled = "#BC3C29FF", cover_labelled = "#E18727FF", unmatch_labelled = "#0072B5FF"
  ))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        plot.title=element_text(color="black", size=10, vjust=0.5),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",linewidth = .4),
        axis.ticks.x=element_line(color="black",linewidth = .4,lineend = 1),
        axis.ticks.y=element_line(color="black",linewidth = .4,lineend = 1),
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 9))+
  labs(x="unilung annotation",y="Percentage of cells",title="Comparison between original and unilung label (level3)")


data = perlv4[-6,]
data = melt(data)
colnames(data) = c('celltype','labelled','percent')

ggplot(data, aes(x = celltype, weight = percent, fill = labelled))+
  geom_bar(position = "stack",width=.4)+
  scale_fill_manual(values=c(
    match_labelled = "#BC3C29FF", cover_labelled = "#E18727FF", unmatch_labelled = "#0072B5FF"
  ))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        plot.title=element_text(color="black", size=10, vjust=0.5),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",linewidth = .4),
        axis.ticks.x=element_line(color="black",linewidth = .4,lineend = 1),
        axis.ticks.y=element_line(color="black",linewidth = .4,lineend = 1),
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 9))+
  labs(x="unilung annotation",y="Percentage of cells",title="Comparison between original and unilung label (level4)")









