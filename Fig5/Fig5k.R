
setwd('/data_hdd4/yuseong/cancer_aging/analysis/h.github/Fig5/')

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

create_custom_vlnplot <- function(df, y_var, y_label, p_y_pos, y_limits = NULL, scientific = TRUE, cols = NULL, comp) {
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



#### malignant
mal.so = qread('data/5k/spotiphy_malignant_seuratobject.qs')

mal.score = read.csv('data/5k/spotiphy_mal_pathway_scores.csv', row.names=1)

mal.so = AddMetaData(mal.so, mal.score)


#### Extended Data Fig 5d - senescence scores in malignant cells
plot_dat = data.frame(group = mal.so$cell_type,
                      SEN.MAYO = mal.so$SEN_MAYO,
                      RESPONSE.TO.UV = mal.so$CELLULAR.RESPONSE.TO.UV,
                      EMT = mal.so$EMT)

plot_dat2 = reshape2::melt(plot_dat, id.vars = 'group', variable.name = 'score')


comp_m <- list(c("Sen_Malignant", "Malignant"))
cols = c('Malignant' = '#6A5ACD', 'Sen_Malignant' = '#E35D5D')

p1 <- create_custom_vlnplot(plot_dat2[plot_dat2$score == 'EMT',], "EMT", "EMT", p_y_pos = 1.0, scientific = TRUE, cols = cols, comp = comp_m)
p2 <- create_custom_vlnplot(plot_dat2[plot_dat2$score == 'RESPONSE.TO.UV',], "RESPONSE.TO.UV", "RESPONSE.TO.UV", p_y_pos = 1.2, scientific = TRUE, cols = cols, comp = comp_m)
p3 <- create_custom_vlnplot(plot_dat2[plot_dat2$score == 'SEN.MAYO',], "SEN.MAYO", "SEN.MAYO", p_y_pos = 0.6, scientific = TRUE, cols = cols, comp = comp_m)

grid.arrange(p1, p2, p3, ncol = 3)


ggsave('figures/FigS5d_mal_senscores_norm.pdf', grid.arrange(p1, p2, p3, ncol = 3),
       width = 10, height=4, dpi=300)


#### Fig 5k top - pathway scores in malignant cells
mal.so2 = subset(mal.so, subset = near3 != 'exclude')

comp_nm <- list(c('near_mal', 'Sen_Malignant'),c('far_mal', 'Sen_Malignant'),
               c('far_mal', 'near_mal'))

mal.so2$near3 = factor(mal.so2$near3, levels = c('Sen_Malignant','near_mal', 'far_mal'))

plot_dat_nm = data.frame(group = mal.so2$near3,
                      EMT = mal.so2$EMT,
                      INVASIVE = mal.so2$INVASIVE,
                      METASTASIS = mal.so2$METASTASIS)

plot_dat_nm2 = reshape2::melt(plot_dat_nm, id.vars = 'group', variable.name = 'score')


cols_nm = c('Sen_Malignant' = 'firebrick2', 'near_mal' = 'deeppink4', 'far_mal' = 'thistle2')

p4 <- create_custom_vlnplot(plot_dat_nm2[plot_dat_nm2$score == 'EMT',], "EMT", "EMT",scientific = TRUE, cols = cols_nm, comp = comp_nm)
p5 <- create_custom_vlnplot(plot_dat_nm2[plot_dat_nm2$score == 'INVASIVE',], "INVASIVE", "INVASIVE", scientific = TRUE, cols = cols_nm, comp = comp_nm)
p6 <- create_custom_vlnplot(plot_dat_nm2[plot_dat_nm2$score == 'METASTASIS',], "METASTASIS", "METASTASIS", scientific = TRUE, cols = cols_nm, comp = comp_nm)

grid.arrange(p4, p5, p6, ncol = 3)

ggsave('figures/Fig5k_st_mal_ssGSEAscore_norm.pdf', grid.arrange(p4, p5, p6, ncol = 3),
       width = 10, height=4, dpi=300)



#### endothelial
end.so = qread('data/5k/spotiphy_endothelial_seuratobject.qs')

end.score = read.csv('data/5k/spotiphy_endo_pathway_scores.csv', row.names=1)

end.so = AddMetaData(end.so, end.score)


#### Extended Data Fig 5e - senescence scores in endothelial cells
plot_dat3 = data.frame(group = end.so$cell_type,
                      SEN.MAYO = end.so$SEN_MAYO,
                      RESPONSE.TO.UV = end.so$CELLULAR.RESPONSE.TO.UV,
                      EndMT = end.so$EndMT)

plot_dat4 = reshape2::melt(plot_dat3, id.vars = 'group', variable.name = 'score')


comp_e <- list(c("Sen_Endothelial", "Endothelial"))
cols_e = c('Endothelial' = '#1F78B4', 'Sen_Endothelial' = '#000000')

p7 <- create_custom_vlnplot(plot_dat4[plot_dat4$score == 'EndMT',], "EndMT", "EndMT", p_y_pos = 1.0, scientific = TRUE, cols = cols_e, comp = comp_e)
p8 <- create_custom_vlnplot(plot_dat4[plot_dat4$score == 'RESPONSE.TO.UV',], "RESPONSE.TO.UV", "RESPONSE.TO.UV", p_y_pos = 1.2, scientific = TRUE, cols = cols_e, comp = comp_e)
p9 <- create_custom_vlnplot(plot_dat4[plot_dat4$score == 'SEN.MAYO',], "SEN.MAYO", "SEN.MAYO", p_y_pos = 0.6, scientific = TRUE, cols = cols_e, comp = comp_e)

grid.arrange(p7, p8, p9, ncol = 3)

ggsave('figures/FigS5e_endo_3scores_norm.pdf', grid.arrange(p1, p2, p3, ncol = 3),
       width = 10, height=4, dpi=300)


#### Fig 5k bottom - pathway scores in endothelial cells
end.so2 = subset(end.so, subset = near3 != 'exclude')

comp_ne <- list(c('near_endo', 'Sen_Endothelial'),c('far_endo', 'Sen_Endothelial'),
               c('far_endo', 'near_endo'))

end.so2$near3 = factor(end.so2$near3, levels = c('Sen_Endothelial','near_endo', 'far_endo'))

plot_dat_ne = data.frame(group = end.so2$near3,
                      PERM = end.so2$PERM,
                      CELL_JUNCTION = end.so2$CELL_JUNCTION,
                      CELL_ADHESION = end.so2$CELL_ADHESION)

plot_dat_ne2 = reshape2::melt(plot_dat_ne, id.vars = 'group', variable.name = 'score')


cols_ne = c('Sen_Endothelial' = 'black', 'near_endo' = 'gray60', 'far_endo' = 'grey88')

p10 <- create_custom_vlnplot(plot_dat_ne2[plot_dat_ne2$score == 'PERM',], "PERM", "PERM", scientific = TRUE, cols = cols_ne, comp = comp_ne)
p11 <- create_custom_vlnplot(plot_dat_ne2[plot_dat_ne2$score == 'CELL_JUNCTION',], "CELL_JUNCTION", "CELL_JUNCTION", scientific = TRUE, cols = cols_ne, comp = comp_ne)
p12 <- create_custom_vlnplot(plot_dat_ne2[plot_dat_ne2$score == 'CELL_ADHESION',], "CELL_ADHESION", "CELL_ADHESION", scientific = TRUE, cols = cols_ne, comp = comp_ne)

grid.arrange(p10, p11, p12, ncol = 3)

ggsave('figures/Fig5k_st_endo_ssGSEAscore_norm.pdf', grid.arrange(p10, p11, p12, ncol = 3),
       width = 10, height=4, dpi=300)

