
setwd('/data_hdd4/yuseong/cancer_aging/analysis/h.github/Fig1')

rm(list=ls())

library(WGCNA)       
library(qs) 
library(dplyr)         
library(tibble)            
library(tidyr)            
library(ggplot2)           
library(ComplexHeatmap)     
library(circlize)           
library(grid)        


### SET CORES TO USE
library(BiocParallel)
library(RhpcBLASctl)
blas_set_num_threads(30)
omp_set_num_threads(30)

enableWGCNAThreads(nThreads = 30)

### PLOT OPTIONS
col_fun2 = colorRamp2(c(-4, 0, 4), c("#313695", "white", "#a50026"))

cell_fun = function(j, i, x, y, width, height, fill) {
  p_val <- p_mat[i, j]
  if (!is.na(p_val) & p_val < 0.05) {
    star_text <- case_when(
      p_val < 0.001 ~ "***",
      p_val < 0.01  ~ "**",
      TRUE          ~ "*"
    )
    grid.text(star_text, x, y, gp = gpar(fontsize = 10, fontface = "bold"))
  }
}


### WGCNA
expr_normalized = qread('data/1f/RTcellline_normalized_expression.qs')

input_mat = t(expr_normalized)
input_mat[1:5,1:10] 

softPower = 14

net = blockwiseModules(input_mat, power = softPower, networkType = 'signed', TOMType = 'signed',
                       minModuleSize = 50, reassignThreshold = 0, 
                       mergeCutHeight = 0.3, corType = "bicor",            
                       maxPOutliers = 0.05,
                       numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE, 
                       saveTOMFileBase = 'data/RTcellline_WGCNA_TOM',
                       nThreads = 30)

qsave(net, 'data/RTcellline_WGCNA_networks_bicor.qs', preset='fast')


module_colors = labels2colors(net$colors)
names(module_colors) = colnames(expression)

module_mapping <- c(
  "MEbrown"     = "M1",
  "MEblack"     = "M2",
  "MEblue"      = "M3",
  "MEgreen"     = "M4",
  "MEyellow"    = "M5",
  "MEred"       = "M6",
  "MEturquoise" = "M7",
  "MEgrey"      = "M8"
)

me = list()
for (col in unique(module_colors)) {
  m.name = module_mapping[paste0('ME', col)]
  me[[m.name]] = names(module_colors)[module_colors ==col]
}

saveRDS(me, 'data/1f/RTcellline_WGCNA_Module_genes_bicor.rds')



### Fig 1f left - difference in GSVA module scores across our celllines
khu.cl = read.csv('data/1f/KHUcellline_WGCNA_GSVA_results.csv', row.names=1)

khu.cl$group = factor(khu.cl$group, levels = c('SEN','NR'))
khu.cl$ME = factor(khu.cl$ME, levels = paste0('M',seq(1:8)))

heatmap_stats <- khu.cl %>%
  group_by(ME) %>%
  summarize(
    t_stat = t.test(Score ~ group)$statistic,
    p_val  = t.test(Score ~ group)$p.value,
    .groups = "drop"
  )

plot_mat <- heatmap_stats %>%
  dplyr::select(ME, t_stat) %>%
  tibble::column_to_rownames("ME") %>%
  as.matrix()

p_mat <- heatmap_stats %>%
  dplyr::select(ME, p_val) %>%
  tibble::column_to_rownames("ME") %>%
  as.matrix()

a = Heatmap(plot_mat, 
            name = "t-statistic", 
            col = col_fun2,
            rect_gp = gpar(col = "white", lwd = 1), 
            column_names_centered = TRUE,
            cluster_rows = F,                    
            cluster_columns = F, 
            column_names_rot = 0,
            cell_fun = cell_fun,                 
            row_names_side = "left",
            row_title = "Gene Modules (GSVA)",
            heatmap_legend_param = list(direction = "horizontal"))
a

# SAVE
pdf('figures/Fig1f_KHUcellline_WGCNA_GSVA_t-statistic.pdf',
    width=3, height=6)

a

dev.off()



### Fig 1f middle - difference in GSVA module scores across independent celllines
rt.cl = read.csv('data/1f/RTcellline_WGCNA_GSVA_results.csv', row.names=1)

rt.cl$group = factor(rt.cl$group, levels = c('RT','NT'))
rt.cl$ME = factor(rt.cl$ME, levels = paste0('M',seq(1:8)))

heatmap_stats2 <- rt.cl %>%
  group_by(ME, cluster) %>%
  summarize(
    t_stat = t.test(Score ~ group)$statistic,
    p_val  = t.test(Score ~ group)$p.value,
    .groups = "drop"
  )

plot_mat2 <- heatmap_stats2 %>%
  dplyr::select(ME, cluster, t_stat) %>%
  tidyr::pivot_wider(names_from = cluster, values_from = t_stat) %>%
  tibble::column_to_rownames("ME") %>%
  as.matrix()

plot_mat2 = plot_mat2[paste0('M', seq(1,8)),]

p_mat2 <- heatmap_stats2 %>%
  dplyr::select(ME, cluster, p_val) %>%
  tidyr::pivot_wider(names_from = cluster, values_from = p_val) %>%
  tibble::column_to_rownames("ME") %>%
  as.matrix()
p_mat2 = p_mat2[paste0('M', seq(1,8)),]


b = Heatmap(plot_mat2, 
            name = "t-statistic", 
            col = col_fun2,
            rect_gp = gpar(col = "white", lwd = 1), 
            column_names_centered = TRUE,
            cluster_rows = F,                    
            cluster_columns = F, 
            column_names_rot = 0,
            cell_fun = cell_fun,                 
            row_names_side = "left",
            row_title = "Gene Modules (GSVA)",
            heatmap_legend_param = list(direction = "horizontal"))
b

# SAVE
pdf('figures/Fig1f_RTcellline_WGCNA_GSVA_t-statistic.pdf',
    width=4, height=6)

b

dev.off()


### Fig 1f right - Pathway enrichment scores in WGCNA module genes
me_pathways = read.table('data/1f/RTcellline_WGCNA_module_genes_pathway_results.tsv', sep='\t',
                         header=T)

target_pathways = c(        
  ## SENESCENCE
  "SAUL_SEN_MAYO",
  "HALLMARK_P53_PATHWAY",
  "GOBP_REGULATION_OF_CELL_DIVISION",
  
  ## SASP
  "GOBP_TYPE_II_INTERFERON_PRODUCTION" ,
  'HALLMARK_INFLAMMATORY_RESPONSE',
  'GOBP_CELLULAR_RESPONSE_TO_INTERLEUKIN_1',
  "GOBP_CYTOKINE_MEDIATED_SIGNALING_PATHWAY" ,
  "REACTOME_INTERLEUKIN_10_SIGNALING",    
  
  ## EMT    
  "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX",  
  "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
  "GOBP_REGULATION_OF_CELL_CELL_ADHESION",  
  "GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION",
  "REACTOME_ECM_PROTEOGLYCANS" ,
  "GOBP_CELL_MATRIX_ADHESION",
  
  ## RESPONSE TO UV
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_HYPOXIA",
  "GOBP_APOPTOTIC_CELL_CLEARANCE"
)

me_pathways$module = factor(me_pathways$module, levels = rev(c('M1','M2','M3','M4','M5','M6','M7','M8'))) 
me_pathways$ID = factor(me_pathways$ID, levels = (target_pathways)) 


ggplot(final_plot_df, aes(x = ID, y = module, fill = foldchange_p)) +
  geom_tile(color = "white", lwd = 0.5) + # 격자선 추가
  scale_fill_gradient(low = "white", high = "#238B45", name = "fold enrichment") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    panel.grid.major = element_blank(),
    legend.position = "right"
  )


ggsave('figures/Fig1f_WGCNA_GENESET_score_bicor.pdf',
       width=9, height = 9, dpi=300)

