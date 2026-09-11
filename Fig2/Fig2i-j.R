
setwd('/data_hdd4/yuseong/cancer_aging/analysis/h.github/Fig2/')

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


### Fig 2i - GSVA scores across malignant clusters
mal.so = qread('data/2e_g_i/patient_mal_seurat_+GSVA+ME_ssgsea_bicor.qs')
mal.cluster = read.csv('data/2e_g_i/patient_mal_clusters.csv', row.names=1)

mal.so = AddMetaData(mal.so, mal.cluster)
mal.so$sen_cluster = factor(mal.so$sen_cluster, levels = c('C0','C1', 'C2', 'C3', 'C4', 'C5','C6', 'C7','C8', 'C9'))

comp_m <- list(c('C0', 'C6'),c('C1', 'C6'),c('C2', 'C6'), c('C3', 'C6'),
               c('C4', 'C6'),c('C5', 'C6'),c('C7', 'C6'),c('C8', 'C6'),c('C9', 'C6'))


plot_dat = data.frame(group = mal.so$sen_cluster,
                      SEN.MAYO = mal.so$SEN_MAYO,
                      RESPONSE.TO.UV = mal.so$CELLULAR.RESPONSE.TO.UV,
                      EMT = mal.so$EMT,
                      M3 = mal.so$M3)

plot_dat2 = reshape2::melt(plot_dat, id.vars = 'group', variable.name = 'score')


cols = c("#A6CEE3", "#1F78B4" ,"#B2DF8A" ,"#33A02C","#FB9A99" ,
         "#FDBF6F","#E31A1C", "#FF7F00", "#CAB2D6","#6A3D9A")

p1 <- create_custom_vlnplot(plot_dat2[plot_dat2$score == 'SEN.MAYO',], "SEN.MAYO", "SEN.MAYO", scientific = TRUE, cols = cols, comp = comp_m)
p2 <- create_custom_vlnplot(plot_dat2[plot_dat2$score == 'RESPONSE.TO.UV',], "RESPONSE.TO.UV", "RESPONSE.TO.UV", scientific = TRUE, cols = cols, comp = comp_m)
p3 <- create_custom_vlnplot(plot_dat2[plot_dat2$score == 'EMT',], "EMT", "EMT",scientific = TRUE, cols = cols, comp = comp_m)
p4 <- create_custom_vlnplot(plot_dat2[plot_dat2$score == 'M3',], "M3", "M3", scientific = TRUE, cols = cols, comp = comp_m)

grid.arrange(p1, p2, p3,p4, ncol = 4)

ggsave('figures/Fig2i_sc_mal_GSVAscore_vlnplt_bicor_norm.pdf', grid.arrange(p1, p2, p3,p4, ncol = 4),
       width = 12, height=4, dpi=300)




### Fig 2j - GSVA scores across endothelial clusters
end.so = qread('data/2f_h_j/patient_endo_seurat_+GSVA+ME_ssgsea_bicor.qs')
end.cl = read.csv('data/2i-j/patient_endo_clusters.csv', row.names=1)

end.so = AddMetaData(end.so, end.cl)
end.so$Endothelial_clusters = factor(end.so$Endothelial_clusters, levels = c('E0','E1', 'E2', 'E3', 'E4', 'E5','E6'))


comp_e <- list(c('E0', 'E5'),c('E1', 'E5'),c('E2', 'E5'), c('E3', 'E5'),
               c('E4', 'E5'))


plot_dat3 = data.frame(group = end.so$Endothelial_clusters,
                      SEN.MAYO = end.so$SEN_MAYO,
                      RESPONSE.TO.UV = end.so$CELLULAR.RESPONSE.TO.UV,
                      EMT = end.so$EMT,
                      M3 = end.so$M3)

plot_dat4 = reshape2::melt(plot_dat3, id.vars = 'group', variable.name = 'score')


cols2 = c("#A6CEE3", "#1F78B4" ,"#B2DF8A" ,"#33A02C","#FB9A99" ,
         "#FDBF6F","#E31A1C", "#FF7F00", "#CAB2D6","#6A3D9A")

p1 <- create_custom_vlnplot(plot_dat4[plot_dat4$score == 'SEN.MAYO',], "SEN.MAYO", "SEN.MAYO", scientific = TRUE, cols = cols2, comp = comp_e)
p2 <- create_custom_vlnplot(plot_dat4[plot_dat4$score == 'RESPONSE.TO.UV',], "RESPONSE.TO.UV", "RESPONSE.TO.UV", scientific = TRUE, cols = cols2, comp = comp_e)
p3 <- create_custom_vlnplot(plot_dat4[plot_dat4$score == 'EMT',], "EMT", "EMT",scientific = TRUE, cols = cols2, comp = comp_e)
p4 <- create_custom_vlnplot(plot_dat4[plot_dat4$score == 'M3',], "M3", "M3", scientific = TRUE, cols = cols2, comp = comp_e)

grid.arrange(p1, p2, p3,p4, ncol = 4)

ggsave('figures/Fig2j_sc_endo_GSVAscore_vlnplt_bicor_norm.pdf', grid.arrange(p1, p2, p3,p4, ncol = 4),
       width = 12, height=4, dpi=300)




