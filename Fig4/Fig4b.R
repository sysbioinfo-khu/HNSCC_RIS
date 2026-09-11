
setwd('/data_hdd4/yuseong/cancer_aging/analysis/h.github/')
rm(list=ls())

library(dplyr)
library(qs)
library(Seurat)
library(SeuratObject)
library(patchwork)
library(escape)
library(BiocParallel)
library(ggpubr)
library(jsonlite)
library(Matrix)
library(gridExtra)

## -------------
## plot function
## -------------

create_custom_vlnplot <- function(df, y_var, y_label, p_y_pos, y_limits = NULL, scientific = TRUE, cols = NULL, comp ) {
  y_axis_label <- bquote(.(y_label) ~ 'score')
  
  p <- ggplot(df, aes(x = group, y = value, color = group)) +
    geom_violin(aes(fill=NULL), trim = FALSE, position = position_dodge(width = 0.5), alpha = 0) +
    
    geom_boxplot(
      width = 0.05, 
      fatten = 0,               
      outlier.shape = NA,       
      aes(fill = group, color = group),       
      size = 0.5             
    ) +
    stat_summary(
      fun = median, 
      geom = "point", 
      shape = 21,            
      size = 1.5, 
      fill = "gray80",        
      color = "black",        
      stroke = 0.5
    ) +
    stat_compare_means(
      comparisons = comp, 
      method = "t.test", 
      label = 'p.signif',
      # label.y = p_y_pos,
      size = 3
    ) +
    
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    
    theme_classic() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(face = "bold", size = 11, color = "black"),
      axis.title.y = element_text(size = 12, color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 13),
      panel.border = element_blank(),
      legend.position = "none"
    ) +
    labs(y = y_axis_label)
  
  if (!is.null(y_limits)) {
    p <- p + scale_y_continuous(limits = y_limits, expand = expansion(mult = c(0.05, 0.15)))
  } else {
    p <- p + scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
  }
  
  return(p)
}

### Extended Data Fig 4b & Fig 4b
mal.so = qread('Fig2/data/2e_g_i/patient_mal_seurat_+GSVA+ME_ssgsea_bicor.qs')
mal.cluster = read.csv('Fig2/data/2e_g_i/patient_mal_clusters.csv', row.names=1)

mal.so = AddMetaData(mal.so, mal.cluster)
mal.so$sen_cluster = factor(mal.so$sen_cluster, levels = c('C0','C1', 'C2', 'C3', 'C4', 'C5','C6', 'C7','C8', 'C9'))

# add TF activities
acts.mat = read.csv('Fig4/data/4b/sc_mal_VIPER_pivot.csv', row.names=1)

mal.so = AddMetaData(mal.so, acts.mat)



comp_m <- list(c('C0', 'C6'),c('C1', 'C6'),c('C2', 'C6'), c('C3', 'C6'),
               c('C4', 'C6'),c('C5', 'C6'),c('C7', 'C6'),c('C8', 'C6'),c('C9', 'C6'))

cols = c("#A6CEE3", "#1F78B4" ,"#B2DF8A" ,"#33A02C","#FB9A99" ,
         "#FDBF6F","#E31A1C", "#FF7F00", "#CAB2D6","#6A3D9A")


### Extended Data Fig 4b - TF expression across malignant clusters
VlnPlot(mal.so, features = c('POU2F2','EGR1', 'FLI1','NFATC2'),pt.size=0, ncol=4, group.by = 'sen_cluster', cols = cols) &
  stat_summary(geom='crossbar', fun=mean, width=0.3, color='black') &
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) &
  stat_compare_means(comparisons = comp_m, method = "t.test", label = "p.signif")

ggsave('Fig4/figures/FigS4b_sc_mal_TFexp_vlnplt.pdf',
       width = 12, height=4, dpi=300)


### Fig 4b - TF activity across malignant clusters
plot_dat = data.frame(group = mal.so$sen_cluster,
                      Regulon_POU2F2 = mal.so$Regulon_POU2F2,
                      Regulon_EGR1 = mal.so$Regulon_EGR1,
                      Regulon_FLI1 = mal.so$Regulon_FLI1, 
                      Regulon_NFATC2 = mal.so$Regulon_NFATC2)

plot_dat2 = reshape2::melt(plot_dat, id.vars = 'group', variable.name = 'score')


p1 <- create_custom_vlnplot(plot_dat2[plot_dat2$score == 'Regulon_POU2F2',], "Regulon_POU2F2", "Regulon_POU2F2", scientific = TRUE, cols = cols, comp = comp_m)
p2 <- create_custom_vlnplot(plot_dat2[plot_dat2$score == 'Regulon_EGR1',], "Regulon_EGR1", "Regulon_EGR1", scientific = TRUE, cols = cols, comp = comp_m)
p3 <- create_custom_vlnplot(plot_dat2[plot_dat2$score == 'Regulon_FLI1',], "Regulon_FLI1", "Regulon_FLI1",scientific = TRUE, cols = cols, comp = comp_m)
p4 <- create_custom_vlnplot(plot_dat2[plot_dat2$score == 'Regulon_NFATC2',], "Regulon_NFATC2", "Regulon_NFATC2", scientific = TRUE, cols = cols, comp = comp_m)


grid.arrange(p1, p2, p3,p4, ncol = 4)

ggsave('Fig4/figures/Fig4b_sc_mal_TFscore_vlnplt.pdf', grid.arrange(p1, p2, p3,p4, ncol = 4),
       width = 12, height=4, dpi=300)




### Extended Data Fig 4g
end.so = qread('Fig2/data/2f_h_j/patient_endo_seurat_+GSVA+ME_ssgsea_bicor.qs')
end.cl = read.csv('Fig2/data/2f_h_j/patient_endo_clusters.csv', row.names=1)

end.so = AddMetaData(end.so, end.cl)
end.so$Endothelial_clusters = factor(end.so$Endothelial_clusters, levels = c('E0','E1', 'E2', 'E3', 'E4', 'E5','E6'))

# add TF activity
acts.mat2 = read.csv('Fig4/data/4b/sc_endo_VIPER_pivot.csv', row.names=1)

end.so = AddMetaData(end.so, acts.mat2)


comp_e <- list(c('E0', 'E5'),c('E1', 'E5'),c('E2', 'E5'), c('E3', 'E5'),
               c('E4', 'E5'))

end_cols = c("#7FC97F", "#BEAED4", "#FDC086", "#FFCC00", "#386CB0", "#000000", "#BF5B17")

### Extended Data Fig 4g top - TF expression across endothelial clusters
VlnPlot(end.so, features = c('POU2F2','EGR1'),pt.size=0, ncol=2, group.by = 'Endothelial_clusters', cols = end_cols) &
  stat_summary(geom='crossbar', fun=mean, width=0.3, color='grey') &
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) &
  stat_compare_means(comparisons = comp_e, method = "t.test", label = "p.signif")

ggsave('Fig4/figures/FigS4g_sc_endo_TFexp_vlnplt.pdf',
       width = 8, height=4, dpi=300)


### Extended Data Fig 4g middle - TF activity across endothelial clusters
plot_dat3 = data.frame(group = end.so$Endothelial_clusters,
                      Regulon_POU2F2 = end.so$Regulon_POU2F2,
                      Regulon_EGR1 = end.so$Regulon_EGR1,
                      Regulon_FLI1 = end.so$Regulon_FLI1,
                      Regulon_NFATC2 = end.so$Regulon_NFATC2)

plot_dat4 = reshape2::melt(plot_dat3, id.vars = 'group', variable.name = 'score')


p5 <- create_custom_vlnplot(plot_dat4[plot_dat4$score == 'Regulon_POU2F2',], "Regulon_POU2F2", "Regulon_POU2F2", scientific = TRUE, cols = end_cols, comp = comp_e)
p6 <- create_custom_vlnplot(plot_dat4[plot_dat4$score == 'Regulon_EGR1',], "Regulon_EGR1", "Regulon_EGR1", scientific = TRUE, cols = end_cols, comp = comp_e)
p7 <- create_custom_vlnplot(plot_dat4[plot_dat4$score == 'Regulon_FLI1',], "Regulon_FLI1", "Regulon_FLI1",scientific = TRUE, cols = end_cols, comp = comp_e)
p8 <- create_custom_vlnplot(plot_dat4[plot_dat4$score == 'Regulon_NFATC2',], "Regulon_NFATC2", "Regulon_NFATC2", scientific = TRUE, cols = end_cols, comp = comp_e)


grid.arrange(p5, p6, p7,p8, ncol = 4)

ggsave('Fig4/figures/FigS4g_sc_endo_TFscore_vlnplt.pdf', grid.arrange(p5, p6, p7,p8, ncol = 4),
       width = 8, height=4, dpi=300)


