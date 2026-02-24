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
library(SPATA2)

#Figure 2F
ST=readRDS('ST.rds')
mycolors=c("#451077FF","#721F81FF","#9F2F7FFF","#CD4071FF","#F1605DFF","#FD9567FF","#FEC98DFF","#FCFDBFFF")#viridis magma
fix.sc <- scale_fill_gradientn(colours= mycolors)
TREM2marker=read.csv('TREM2marker.csv')

TREM2marker=list(TREM2marker)
names(TREM2marker)=TREM2marker
ST=AddModuleScore(ST,features = TREM2marker,name = 'TREM2mac')

celltype=c("TREM2mac1")

for (i in celltype) {
  p1 <- SpatialPlot(ST, features = i, combine=FALSE, crop=F, pt.size=1.3)
  p2 <- lapply(p1, function(x) x + fix.sc)
  p3 <- aplot::plot_list(gglist = p2, ncol = 2) + theme(legend.position = "right")
  ggsave(p3, file = paste0('ST',i, '.pdf'), height = 5.5, width = 5.5)
}


spata_obj=readRDS('spata_obj.rds')
p <- plotStsLineplot(
  object = spata_obj,
  id = "T1",
  variables = "TREM2mac1",
  smooth_span = 0.4,
  smooth_se = TRUE,
  line_size = 1.5
)

ggsave(plot=p,file='SPATA_TREM2mac_T1.pdf',height = 3.5,width=3.5)


p1=plotSurface(spata_obj, color_by ="TREM2mac1", pt_size =1.5)
ggsave(plot=p1,file='SPATA_TREM2mac1.pdf',height = 5.5,width=5.5)

#Figure 2G
enrich_result <- read.csv('TREM2TAM_enrich.csv')
enrich_result$pathway <- factor(enrich_result$Description, levels = rev(enrich_result$Description))
p <- ggplot(data = enrich_result,
            aes(x = -log10(pvalue), y = pathway, fill = Count)) +
  geom_bar(stat = "identity",width = 0.8) +
  scale_fill_distiller(palette = "Oranges",direction = 1)  +
  labs(x = "-log10(pvalue)",
       y = "pathway",
       title = "enrichment barplot") +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11))
p
ggsave(plot = p,filename = "TREM2TAM_enrich.pdf",height=2.5,width=6)



