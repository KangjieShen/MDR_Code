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
library(scales)
library(survival)
library(survminer)
library(broom)
library(forestplot)

#Figure 8C
df =read.csv("pancancer_cor.csv")

df_plot <- df %>%
  mutate(
    neglog10P = -log10(P),
    cell = factor(cell, levels = c("AXL+tumor", "CXCL14+CAF", "CXCL13+CD8T")),
    cancer = factor(cancer, levels = c("MM","HCC","NSCLC","ccRCC"))
  )

p <- ggplot(df_plot, aes(x = cell, y = cancer)) +
  geom_point(aes(size = neglog10P, fill = R), shape = 21, color = "white", stroke = 0.25) +
  scale_size_continuous(
    name = expression(-log[10](P)),
    range = c(2, 8)
  ) +
  scale_fill_gradient2(
    name = "R",
    low = "#08519C", mid = "white", high = "#99000D",
    limits = c(-1, 1), oob = squish
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "right"
  )

print(p)

p1=p + geom_text(aes(label = sprintf("%.2f", R)), size = 3, vjust = 0.5)
p1
ggsave(plot=p1,file='pancancer_cor.pdf',height=4,width = 3.5)


#Figure 8D
rt = read.table('clinical_input.txt', sep='\t', header=T, row.names=1)

cor_matrix <- cor(rt[, c("CD68TREM2", "CXCL14CAF", "MELANAAXL", "CD8CXCL13")])
print(cor_matrix)
uniTab = data.frame()

for(i in colnames(rt)[3:ncol(rt)]){ 
  cox <- coxph(as.formula(paste0("Surv(futime, fustat) ~ `", i, "`")), data = rt)
  coxSummary = summary(cox)
  uniTab = rbind(uniTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
uniTab <- as.data.frame(uniTab)

hrtable <- rbind(c("UniCox", NA, NA, NA, NA),
                 uniTab)

tabletext <- cbind(c("Variable", hrtable$id),
                   c("HR", format(round(as.numeric(hrtable$HR), 3), nsmall = 3)),
                   c("HR.95L", format(round(as.numeric(hrtable$HR.95L), 3), nsmall = 3)),
                   c("HR.95H", format(round(as.numeric(hrtable$HR.95H), 3), nsmall = 3)),
                   c("pvalue", formatC(as.numeric(hrtable$pvalue), format = "e", digits = 2)))

tabletext[2,] <- c("UniCox", NA, NA, NA, NA) 


bottom_line <- as.character(nrow(tabletext) + 1)
lines_list <- list("1" = gpar(lwd=2, col="black", lty=1),        
                   "2" = gpar(lwd=1, col="grey50", lty=2))       
lines_list[[bottom_line]] <- gpar(lwd=2, col="black")           


pdf("indep_uniCox.pdf", width = 8, height = max(6, nrow(tabletext)*0.45)) 
forestplot(labeltext = tabletext,
           mean = c(NA, log2(as.numeric(hrtable$HR))),
           lower = c(NA, log2(as.numeric(hrtable$HR.95L))),
           upper = c(NA, log2(as.numeric(hrtable$HR.95H))),
           graph.pos = 6,
           graphwidth = unit(.25, "npc"),
           fn.ci_norm = "fpDrawDiamondCI",
           col = fpColors(box = "#E64B35FF", lines = "#E64B35FF", zero = "black"),
           boxsize = 0.4,
           lwd.ci = 1,
           ci.vertices.height = 0.1, ci.vertices = F,
           zero = 0,
           lwd.zero = 2,
           xticks = c(-1, 0, 1, 2, 3, 4, 5),
           lwd.xaxis = 2,
           xlab = expression("log"[2]~"HR"),
           hrzl_lines = lines_list, 
           txt_gp = fpTxtGp(label = gpar(cex=1.2),
                            ticks = gpar(cex=0.85),
                            xlab = gpar(cex=1),
                            title = gpar(cex=1.5)),
           lineheight = unit(.75, "cm"),
           colgap = unit(0.3, "cm"),
           mar = unit(rep(1.5, times = 4), "cm"),
           new_page = F
)
invisible(dev.off())
