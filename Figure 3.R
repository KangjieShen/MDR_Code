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
library(CellChat)
library(pheatmap)


#Figure 3A
proportions=read.csv('proportions.csv')
cor_result <- cor.test(proportions$TREM2TAM, proportions$CXCL13CD8, method = "spearman")
print(cor_result)

# 绘制并保存相关性散点图
cor_plot <- ggplot(proportions, aes(x = CXCL13CD8, y = TREM2TAM)) +
  geom_point(alpha = 0.5, color = "#fdbb84") +
  geom_smooth(method = "lm", color = "black") +
  theme_classic() +
  labs(x = "CXCL13CD8", y = "TREM2TAM", 
       title = paste0("R = ", round(cor_result$estimate, 2), 
                      ", p = ", formatC(cor_result$p.value, format = "e", digits = 2)))

ggsave("TREM2TAM_CXCL13CD8_Correlation.pdf", plot = cor_plot, width = 5, height = 5)

#Figure 3F
cellchat=readRDS("cellchat.rds")
mycolors=c('white',"#FFF7EC","#FEE8C8","#FDD49E","#FDBB84","#FC8D59","#EF6548","#D7301F","#B30000","#7F0000")

pheatmap::pheatmap(cellchat@net$count, 
                   border_color = "black", 
                   cluster_cols = F, 
                   fontsize = 10, 
                   cluster_rows = F,
                   display_numbers = T,
                   number_color = "black",
                   number_format = "%.0f",
                   color = mycolors,
                   filename = "heatmap_count.pdf",
                   height = 6, 
                   width = 6)


#Figure 3G
heatmap = read.csv('pathways.csv', row.names = 1)
heatmap_log = log10(heatmap) 
mycolors = c("#000004FF","#180F3EFF","#451077FF","#721F81FF","#9F2F7FFF",
              "#CD4071FF","#F1605DFF","#FD9567FF","#FEC98DFF","#FCFDBFFF")

pheatmap(heatmap_log,
         color = colorRampPalette(mycolors)(50),
         cluster_cols = FALSE,
         cluster_rows = FALSE, 
         border_color = 'white',
         cellwidth = 12,
         cellheight = 12,
         angle_col = '90',
         scale = 'row',
         filename = "heatmap.pdf",
         height = 7,
         width = 4)

#Figure 3H 
pathways.show <- c("CD86") 

pdf("Figure3H.pdf", width = 6, height = 6)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord",sources.use = 14, targets.use = c(1:13,15:19))
dev.off()

#Figure 3I
pdf("Figure3I.pdf", width = 5, height = 4)
netAnalysis_contribution(cellchat, signaling = 'CD86')
dev.off()

#Figure 3J
pairLR <- extractEnrichedLR(cellchat, signaling = 'CD86', geneLR.return = FALSE)
pdf("Figure3J.pdf", width = 5, height = 4)
netVisual_individual(cellchat, signaling = 'CD86', pairLR.use = "CD86_CTLA4", layout = "circle",
                     sources.use = 14, targets.use = c(1:13,15:19))
dev.off()
