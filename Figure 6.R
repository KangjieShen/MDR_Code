rm(list=ls())
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(stringr)  
library(ggplot2)
library(RColorBrewer) 
library(ggrepel)
library(aplot)
library(pheatmap)
library(SPATA2)



#Figure 6A
CAF=readRDS('CAF.rds')
mycolors=c("#FEE8C8","#FDD49E","#FDBB84","#FC8D59","#EF6548","#D7301F","#B30000","#7F0000")

heatmap_gene=c('MME','PDPN','TMEM158','NDRG1',
               'CFD','C3','CXCL14','CXCL12',
               'POSTN','MMP11','COL11A1','MMP14')
heatmap_AveE <- AverageExpression(CAF, assays = "SCT", features = heatmap_gene,verbose = TRUE) %>% .$SCT
heatmap_AveE=t(heatmap_AveE)

pheatmap(heatmap_AveE,
         color = colorRampPalette(mycolors)(50),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         border_color = 'white',
         cellwidth = 10,
         cellheight = 10,
         angle_col = '90',
         scale = 'column',
         filename = "htmap_genesets_CAF.pdf",
         height = 3,
         width = 11.5)

#Figure 6D
enrich_result <- read.csv('CXCL14CAF_enrich.csv')
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
ggsave(plot = p,filename = "CXCL14CAF_enrich.pdf",height=2.5,width=6)

#Figure 6E
proportions=read.csv('proportions.csv')
cor_result <- cor.test(proportions$TREM2TAM, proportions$CXCL14CAF, method = "spearman")
print(cor_result)

# 绘制并保存相关性散点图
cor_plot <- ggplot(proportions, aes(x = CXCL14CAF, y = TREM2TAM)) +
  geom_point(alpha = 0.5, color = "#fdbb84") +
  geom_smooth(method = "lm", color = "#fdbb84") +
  theme_classic() +
  labs(x = "CXCL14CAF", y = "TREM2TAM", 
       title = paste0("R = ", round(cor_result$estimate, 2), 
                      ", p = ", formatC(cor_result$p.value, format = "e", digits = 2)))

ggsave("TREM2TAM_CXCL14CAF_Correlation.pdf", plot = cor_plot, width = 5, height = 5)


#Figure 6F
spata_obj=readRDS('spata_obj.rds')
p <- plotStsLineplot(
  object = spata_obj,
  id = "T1",
  variables = c("CXCL14CAF1","TREM2mac1"),
  smooth_span = 0.4,
  smooth_se = TRUE,
  line_size = 1.5
)

ggsave(plot=p,file='SPATA_CXCL14CAF_TREM2mac1_T1.pdf',height = 3.5,width=3.5)


p1=plotSurface(spata_obj, color_by ="CXCL14CAF1", pt_size =1.5)
ggsave(plot=p1,file='SPATA_CXCL14CAF1.pdf',height = 5.5,width=5.5)

