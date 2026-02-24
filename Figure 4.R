rm(list=ls())
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(stringr)  
library(ggplot2)
library(pheatmap)
library(SPATA2)
library(RColorBrewer) 

#Figure 4A
tumor=readRDS('tumor.rds')
Idents(tumor)='tumorcelltype'
cluster_cols=c('Melanogenesis'='#CCEA9E','Immunoregulation'='#ABD0E6',
               'ECM'='#E8C56B','Dedifferentiated'='#E63A46','Neural cell like'='#E19E37')

plot1=DimPlot(tumor, reduction = "umap", group.by = "tumorcelltype",label = T,cols = cluster_cols)
plot1
ggsave(plot1,file='tumor_celltype_umap.pdf',height=4.5,width=6)



#Figure 4B
mycolors=c("#000004FF","#180F3EFF","#451077FF","#721F81FF","#9F2F7FFF","#CD4071FF","#F1605DFF","#FD9567FF","#FEC98DFF","#FCFDBFFF")#viridis magma
heatmap_gene=c('MITF',
               'SOX10',
               'PMEL',
               'TYRP1',
               'HLA-J',
               'HLA-F',
               'GATA3',
               'CCL5',
               'MMP28',
               'MMP20',
               'COL25A1',
               'COL7A1',
               'AXL',
               'ALDH2',
               'CD34',
               'AQP1',
               'NGFR',
               'EGFR',
               'ERBB3'
)

heatmap_AveE <- AverageExpression(tumor, assays = "SCT", features = heatmap_gene,verbose = TRUE) %>% .$SCT
heatmap_AveE=as.matrix(heatmap_AveE)
heatmap_AveE=t(heatmap_AveE)

pdf("markergene_tumor_celltype.pdf")
pheatmap(heatmap_AveE,color=colorRampPalette(mycolors)(50),cluster_cols=F,cluster_rows = F,
         border_color='white',cellwidth=15,cellheight=15,angle_col='45',scale = 'column')
dev.off()

#Figure 4D
spata_obj=readRDS('spata_obj.rds')
p <- plotStsLineplot(
  object = spata_obj,
  id = "T1",
  variables = "Dedifferentiated1",
  smooth_span = 0.4,
  smooth_se = TRUE,
  line_size = 1.5
)

ggsave(plot=p,file='SPATA_Dedifferentiated1_T1.pdf',height = 3.5,width=3.5)


p1=plotSurface(spata_obj, color_by ="Dedifferentiated1", pt_size =1.5)
ggsave(plot=p1,file='SPATA_Dedifferentiated1.pdf',height = 5.5,width=5.5)



#Figure 4E
enrich_result <- read.csv('Dedifferentiated_enrich.csv')
enrich_result$pathway <- factor(enrich_result$Description, levels = rev(enrich_result$Description))
p <- ggplot(data = enrich_result,
            aes(x = -log10(pvalue), y = pathway, fill = Count)) +
  geom_bar(stat = "identity",width = 0.8) +
  scale_fill_distiller(palette = "OrRd",direction = 1)  +
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
ggsave(plot = p,filename = "Dedifferentiated_enrich.pdf",height=2.5,width=6)


#Figure 4G
ST_tumor=readRDS('ST_tumor.rds')
mycolors=c("#451077FF","#721F81FF","#9F2F7FFF","#CD4071FF","#F1605DFF","#FD9567FF","#FEC98DFF","#FCFDBFFF")#viridis magma
fix.sc <- scale_fill_gradientn(colours= mycolors)


p1 <- SpatialPlot(ST_tumor, features = 'AXL', combine=FALSE, crop=F, pt.size=1.3)
p2 <- lapply(p1, function(x) x + fix.sc)
p3 <- aplot::plot_list(gglist = p2, ncol = 2) + theme(legend.position = "right")
p3
ggsave(p3, file = 'ST_tumor_AXL.pdf', height = 5.5, width = 5.5)



