rm(list=ls())
library(Seurat)
library(ggplot2)
library(dplyr)
library(RANN)
library(pheatmap)
library(readr)
library(ggpubr)


xenium.obj=readRDS('MDR_MM.rds')

#Figure 7B
SpatialSmoothing <- function(coords, expr_mat, gene_name, k = 20) {
  neighbors <- nn2(data = coords, query = coords, k = k)
  gene_vals <- expr_mat[gene_name, ]
  smoothed_vals <- apply(neighbors$nn.idx, 1, function(idx) mean(gene_vals[idx]))
  return(smoothed_vals)
}

DefineSpatialNexus <- function(obj, axl_cutoff = 0.1, trem2_cutoff = 0.1, cxcl12_cutoff = 0.1, dist_cutoff = 50) {
  
  
  coords_df <- GetTissueCoordinates(obj)
  
  
  if ("cell" %in% colnames(coords_df)) {
    rownames(coords_df) <- coords_df$cell
  } else if (ncol(coords_df) >= 3) {
    rownames(coords_df) <- coords_df[, 3]
  }
  
  coords <- coords_df[, 1:2]
  expr_data <- GetAssayData(obj, layer = "data")
  
  obj$AXL_smooth <- SpatialSmoothing(coords, expr_data, "AXL", k = 20)
  obj$TREM2_smooth <- SpatialSmoothing(coords, expr_data, "TREM2", k = 20)
  obj$CXCL12_smooth <- SpatialSmoothing(coords, expr_data, "CXCL12", k = 20)
  
  ids_tumor <- WhichCells(obj, expression = major_cell_type == "Tumor")
  ids_tam <- WhichCells(obj, expression = major_cell_type == "TAM")
  
  dist_res <- nn2(data = coords[ids_tam, ], query = coords[ids_tumor, ], k = 1)
  
  obj$dist_to_tam <- Inf
  obj$dist_to_tam[ids_tumor] <- dist_res$nn.dists
  
  is_nexus_tumor <- (obj$major_cell_type == "Tumor") & (obj$AXL_smooth > axl_cutoff) & (obj$dist_to_tam < dist_cutoff)
  is_nexus_tam <- (obj$major_cell_type == "TAM") & (obj$TREM2_smooth > trem2_cutoff)
  is_nexus_caf <- (obj$major_cell_type == "CAF") & (obj$CXCL12_smooth > cxcl12_cutoff)
  
  obj$nexus_identity <- "Other"
  obj$nexus_identity[obj$major_cell_type == "T cell"] <- "T cell"
  obj$nexus_identity[obj$major_cell_type == "Tumor"] <- "AXL- Tumor"
  
  obj$nexus_identity[is_nexus_tumor] <- "AXL+ Tumor"
  obj$nexus_identity[is_nexus_tam]   <- "TREM2+ TAM"
  obj$nexus_identity[is_nexus_caf]   <- "CXCL12+ CAF"
  
  return(obj)
}


xenium.obj <- DefineSpatialNexus(xenium.obj)


xenium.obj$nexus_identity <- factor(
  xenium.obj$nexus_identity, 
  levels = c("AXL- Tumor", "AXL+ Tumor", "T cell", "Other", "TREM2+ TAM", "CXCL12+ CAF")
)


nexus_colors_final <- c(
  "AXL- Tumor"    = "#f76f76",  
  "AXL+ Tumor"    = "#6fad6f",  
  "T cell"        = "#CCEA9E",  
  "Other"         = "#EFEFEF",  
  "TREM2+ TAM"   = "#FF00FF",  
  "CXCL12+ CAF"  = "#ff8000"   
)


p_final <- ImageDimPlot(xenium.obj, 
                        group.by = "nexus_identity", 
                        cols = nexus_colors_final,
                        dark.background = FALSE, 
                        boundaries = "segmentations",
                        border.size = NA,
                        axes = FALSE)


p1_fixed <- p_final + 
  guides(fill = guide_legend(title = "nexus_identity", 
                             override.aes = list(size = 5, alpha = 1, linewidth = 1)))
p1_fixed
ggsave(filename = 'Spatial_celltype_clean.pdf', plot = p1_fixed, height = 5.5, width = 5.5)

#Figure 7C
IdentifyCellularNeighborhoods <- function(obj, 
                                          cell_type_col = "major_cell_type", 
                                          k_neighbors = 20, 
                                          k_clusters = 6) {
  
  coords_df <- GetTissueCoordinates(obj)
  if ("cell" %in% colnames(coords_df)) {
    rownames(coords_df) <- coords_df$cell
  } else if (ncol(coords_df) >= 3) {
    rownames(coords_df) <- coords_df[, 3]
  }
  coords <- coords_df[, 1:2]
  
  cell_types <- obj[[cell_type_col]][, 1]
  type_levels <- sort(unique(cell_types))
  
  knn.res <- nn2(data = coords, query = coords, k = k_neighbors + 1)
  neighbors <- knn.res$nn.idx[, 2:(k_neighbors + 1)]
  
  cn_matrix <- matrix(0, nrow = nrow(coords), ncol = length(type_levels))
  colnames(cn_matrix) <- type_levels
  
  target_types <- cell_types[neighbors]
  dim(target_types) <- dim(neighbors)
  
  for(i in 1:nrow(coords)) {
    tab <- table(target_types[i, ])
    cn_matrix[i, names(tab)] <- as.numeric(tab) / k_neighbors
  }
  
  set.seed(123) 
  km_res <- kmeans(cn_matrix, centers = k_clusters, iter.max = 100, nstart = 25)
  
  obj$CN_ID <- paste0("CN_", km_res$cluster)
  
  return(obj)
}

xenium.obj <- IdentifyCellularNeighborhoods(xenium.obj, 
                                            cell_type_col = "major_cell_type", 
                                            k_neighbors = 20, 
                                            k_clusters = 6)

cn_palette <- c(
  "CN_1" = "#6fad6f",  
  "CN_4" = "#f76f76",  
  "CN_2" = "#ff8000",  
  "CN_5" = "#3589C1", 
  "CN_3" = "#ABD0E6", 
  "CN_6" = "#CCEA9E" 
)


p_cn <- ImageDimPlot(xenium.obj,
                     group.by = "CN_ID",
                     cols = cn_palette,
                     dark.background = FALSE,
                     boundaries = "segmentations",
                     border.size = NA,
                     axes = FALSE) +
  ggtitle("Unsupervised Cellular Neighborhoods Map")

p_cn_fixed <- p_cn +
  guides(fill = guide_legend(title = "Neighborhood", 
                             override.aes = list(size = 5, alpha = 1, linewidth = 1)))
p_cn_fixed
ggsave(plot = p_cn_fixed, filename = 'CN_Spatial_distribution.pdf', height = 5.5, width = 5.5)

#Figure 7D
CalculateNicheResidence <- function(obj, 
                                    markers = c("AXL", "TREM2", "CXCL12"), 
                                    count_layer = "counts") {
  
  raw_counts <- FetchData(obj, vars = markers, layer = count_layer)

  obj$stat_identity <- "Other"
  
  is_tam <- obj$major_cell_type == "TAM"
  obj$stat_identity[is_tam] <- "TREM2- TAM"
  obj$stat_identity[is_tam & raw_counts$TREM2 > 0] <- "TREM2+ TAM"
  
  is_caf <- obj$major_cell_type == "CAF"
  obj$stat_identity[is_caf] <- "CXCL12- CAF"
  obj$stat_identity[is_caf & raw_counts$CXCL12 > 0] <- "CXCL12+ CAF"
  
  is_tumor <- obj$major_cell_type == "Tumor"
  obj$stat_identity[is_tumor] <- "AXL- Tumor"
  obj$stat_identity[is_tumor & raw_counts$AXL > 0] <- "AXL+ Tumor"
  
  target_cells <- c("AXL+ Tumor", "AXL- Tumor", "TREM2+ TAM", "TREM2- TAM", "CXCL12+ CAF", "CXCL12- CAF")
  
  df_residence <- FetchData(obj, vars = c("CN_ID", "stat_identity")) %>%
    filter(stat_identity %in% target_cells) %>%
    group_by(stat_identity, CN_ID) %>%    
    summarise(Count = n(), .groups = 'drop') %>%
    group_by(stat_identity) %>%
    mutate(Freq = Count / sum(Count)) %>% 
    ungroup()
  
  df_residence$stat_identity <- factor(df_residence$stat_identity, levels = target_cells)
  
  return(df_residence)
}

df_residence <- CalculateNicheResidence(xenium.obj, count_layer = "counts")

p_residence <- ggplot(df_residence, aes(x = stat_identity, y = Freq, fill = CN_ID)) +
  geom_bar(stat = "identity", position = "fill", width = 0.75, color = "black", size = 0.2) +
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  scale_fill_manual(values = cn_palette) +
  theme_classic() +
  labs(y = "Fraction of Cells Residing in Niche", x = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold", color = "black"))

ggsave("CN_Residence_Plot.pdf", p_residence, width = 6, height = 5)

#Figure 7H
sample_info <- data.frame(
  FileName = c("Naive_MM_data.csv", "MDR_MM_data.csv", "MDR_HCC_data.csv"),
  CancerType = c("Naive MM", "MDR MM", "MDR HCC")
)

all_data <- bind_rows(lapply(1:nrow(sample_info), function(i) {
  temp <- read_csv(sample_info$FileName[i], show_col_types = FALSE)
  temp$CancerType <- sample_info$CancerType[i]
  return(temp)
}))

all_data$CancerType <- factor(all_data$CancerType, levels = rev(c("Naive MM", "MDR MM", "MDR HCC")))

p_val_data <- all_data %>%
  group_by(CancerType, Compartment) %>%
  summarise(
    p.val = wilcox.test(Distance ~ Identity)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p.signif = case_when(
      p.val < 0.0001 ~ "****",
      p.val < 0.001 ~ "***",
      p.val < 0.01 ~ "**",
      p.val < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

plot_data_summary <- all_data %>%
  group_by(CancerType, Compartment, Identity) %>%
  summarise(Median_Dist = median(Distance, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Identity, values_from = Median_Dist) %>%
  left_join(p_val_data, by = c("CancerType", "Compartment"))

colors_palette <- c(
  "AXL- Tumor"    = "#f76f76", 
  "AXL+ Tumor"    = "#6fad6f", 
  "CXCL12- CAF"   = "#ABD0E6", 
  "CXCL12+ CAF"   = "#ff8000"  
)

p_dumb_stat <- ggplot(plot_data_summary) +
  
  geom_segment(data = filter(plot_data_summary, Compartment == "Tumor Compartment"),
               aes(x = `AXL- Tumor`, xend = `AXL+ Tumor`, y = CancerType, yend = CancerType),
               color = "#E0E0E0", size = 1.5) +
  
  geom_segment(data = filter(plot_data_summary, Compartment == "Stroma Compartment"),
               aes(x = `CXCL12- CAF`, xend = `CXCL12+ CAF`, y = CancerType, yend = CancerType),
               color = "#E0E0E0", size = 1.5) +
  
  geom_point(aes(x = `AXL- Tumor`, y = CancerType), color = colors_palette["AXL- Tumor"], size = 5, na.rm = TRUE) +
  geom_point(aes(x = `CXCL12- CAF`, y = CancerType), color = colors_palette["CXCL12- CAF"], size = 5, na.rm = TRUE) +
  geom_point(aes(x = `AXL+ Tumor`, y = CancerType), color = colors_palette["AXL+ Tumor"], size = 5, na.rm = TRUE) +
  geom_point(aes(x = `CXCL12+ CAF`, y = CancerType), color = colors_palette["CXCL12+ CAF"], size = 5, na.rm = TRUE) +
  
  geom_text(data = filter(plot_data_summary, Compartment == "Tumor Compartment"),
            aes(x = (`AXL- Tumor` + `AXL+ Tumor`) / 2, y = CancerType, label = p.signif),
            vjust = -0.8, size = 4, fontface = "bold", color = "black") +
  
  geom_text(data = filter(plot_data_summary, Compartment == "Stroma Compartment"),
            aes(x = (`CXCL12- CAF` + `CXCL12+ CAF`) / 2, y = CancerType, label = p.signif),
            vjust = -0.8, size = 4, fontface = "bold", color = "black") +
  
  facet_wrap(~Compartment, scales = "free_x") +
  theme_minimal(base_size = 14) +
  labs(x = "Median Distance to Nearest TREM2+ TAM (µm)", y = "") +
  theme(
    axis.text.y = element_text(face = "bold", color = "black"),
    axis.text.x = element_text(color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 13),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.ticks.x = element_line(color = "black")
  )

print(p_dumb_stat)
ggsave("Figure_Dumbbell_With_Stats.pdf", p_dumb_stat, width = 9, height = 3.5)


