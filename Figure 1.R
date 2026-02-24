rm(list=ls())
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(stringr)  
library(ggplot2)
library(reshape2) 
library(scales)   
library(aplot)

#Figure 1C
all=readRDS('all.rds')
Idents(all)='major_celltype'
celltype1col=c('B'='#88D0A5','CAF'='#ABD0E6','Endothelial'='#3589C1','Epithelial'='#F5FBB2','Mast'='#FEE99B','Melanoma'='#FDBA6A',
               'Myeloid'='#CCEA9E','Plasma'='#4D62AC','SMC'='#68ADD6','TNK'='#F2A361')

plot1=DimPlot(all, reduction = "umap", group.by = "major_celltype",cols = celltype1col,label = T)
plot1
ggsave(plot=plot1,'major_celltype_umap.pdf',height=10,width=11)

plot1=DimPlot(all, reduction = "umap", group.by = "major_celltype",cols = celltype1col,label = T,split.by = "group")
plot1
ggsave(plot=plot1,'major_celltype_umap_split.pdf',height=10,width=20)

pB2_df <- table(all@meta.data$group,all@meta.data$major_celltype) %>% melt()
colnames(pB2_df) <- c("Cluster","Sample","Number")
sample_color <- celltype1col
pB2 <- ggplot(data = pB2_df, aes(x = Cluster, y = Number, fill = Sample)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=sample_color) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  scale_y_continuous(position = "right",labels = percent)+  
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
pB2
ggsave('prop_major_celltype_group_bar.pdf', width = 3, height=4)


#Figure 1D
ST=readRDS('ST.rds')
mycolors=c("#451077FF","#721F81FF","#9F2F7FFF","#CD4071FF","#F1605DFF","#FD9567FF","#FEC98DFF","#FCFDBFFF")#viridis magma
markers=read.csv('major_marker_for_scoring.csv')
fix.sc <- scale_fill_gradientn(colours= mycolors)

CAF=list(markers$CAF)
names(CAF)=CAF
ST=AddModuleScore(ST,features = CAF,name = 'CAF_score')

Myeloid=list(markers$Myeloid)
names(Myeloid)=Myeloid
ST=AddModuleScore(ST,features = Myeloid,name = 'Myeloid_score')

Melanoma=list(markers$Melanoma)
names(Melanoma)=Melanoma
ST=AddModuleScore(ST,features = Melanoma,name = 'Melanoma_score')

celltype=c("CAF_score1","Myeloid_score1",'Melanoma_score1')

for (i in celltype) {
  p1 <- SpatialPlot(ST, features = i, combine=FALSE, crop=F, pt.size=1.3)
  p2 <- lapply(p1, function(x) x + fix.sc)
  p3 <- aplot::plot_list(gglist = p2, ncol = 2) + theme(legend.position = "right")
  ggsave(p3, file = paste0(i, '.pdf'), height = 5.5, width = 5.5)
}

cor_result <- cor.test(ST$CAF_score1, ST$Myeloid_score1, method = "spearman")
print(cor_result)

# 绘制并保存相关性散点图
cor_plot <- ggplot(ST@meta.data, aes(x = Myeloid_score1, y = CAF_score1)) +
  geom_point(alpha = 0.5, color = "#ccea9e") +
  geom_smooth(method = "lm", color = "black") +
  theme_classic() +
  labs(x = "Myeloid", y = "CAF", 
       title = paste0("R = ", round(cor_result$estimate, 2), 
                      ", p = ", formatC(cor_result$p.value, format = "e", digits = 2)))

ggsave("CAF_Myeloid_Correlation.pdf", plot = cor_plot, width = 5, height = 5)