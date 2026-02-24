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

#Figure 5E
mac=readRDS('mac.rds')
df=FindMarkers(mac, only.pos = F,ident.1 ='Mac_TREM2',ident.2 = 'Mac_SPP1',min.pct = 0, logfc.threshold = 0)
df$gene=rownames(df)
df$p_val <-  -log10(df$p_val)
label_df = read.csv('label.csv')
target_genes = label_df[,1] 
df$label_text = ifelse(df$gene %in% target_genes, df$gene, "")

p1 = ggplot(df, aes(avg_log2FC, p_val))+
  geom_hline(yintercept = 0, linetype="dashed", color = "#999999")+  
  geom_vline(xintercept= c(-0.1,0.1), linetype="dashed", color = "#999999")+
  geom_point(aes(size=p_val, color= p_val))+
  scale_color_gradientn(values = seq(0,1,0.2),
                        colors = c("#39489f","#39bbec",'white',"#f9ed36","#f38466","#b81f25"))+
  scale_size_continuous(range = c(1,4))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = c(0.01,0.7),
        legend.justification = c(0,1))+
  guides(col = guide_colourbar(title = "-Log10_P-value"),
         size = "none")+

  geom_text_repel(aes(label=label_text, color = p_val), size = 3, vjust = -5, hjust=5, max.overlaps = 100)+
  xlab("log2 (fold change)")+
  ylab("-log10 (P_value)")

ggsave(filename = 'volcano.pdf', plot = p1, height = 5.5, width = 5.5)


#Figure 5F
data=read.csv('data.csv')
cor_result <- cor.test(data$GAS6, data$TREM2, method = "spearman")
print(cor_result)

cor_plot <- ggplot(data, aes(x = TREM2, y = GAS6)) +
  geom_point(alpha = 0.5, color = "#fdbb84") +
  geom_smooth(method = "lm", color = "black") +
  theme_classic() +
  labs(x = "TREM2", y = "GAS6", 
       title = paste0("R = ", round(cor_result$estimate, 2), 
                      ", p = ", formatC(cor_result$p.value, format = "e", digits = 2)))

ggsave("TREM2_GAS6_Correlation.pdf", plot = cor_plot, width = 5, height = 5)

#Figure 5G
ST=readRDS('ST.rds')
mycolors=c("#451077FF","#721F81FF","#9F2F7FFF","#CD4071FF","#F1605DFF","#FD9567FF","#FEC98DFF","#FCFDBFFF")#viridis magma
fix.sc <- scale_fill_gradientn(colours= mycolors)

p1 <- SpatialPlot(ST, features = 'GAS6', combine=FALSE, crop=F, pt.size=1.3)
p2 <- lapply(p1, function(x) x + fix.sc)
p3 <- aplot::plot_list(gglist = p2, ncol = 2) + theme(legend.position = "right")
p3
ggsave(p3, file = 'ST_GAS6.pdf', height = 5.5, width = 5.5)
