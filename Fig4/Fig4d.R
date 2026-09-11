setwd('/data_hdd4/yuseong/cancer_aging/analysis/h.github/Fig4')

rm(list = ls())


library(dplyr)
library(patchwork)
library(stringr)
library(ggplot2)
library(reshape2)
library(rstatix)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)


### Fig 4d - TF activity in our celllines
tf_acts = read.csv('data/4d/KHUcellline_TFactivities.csv')
reg = c('POU2F2','EGR1','FLI1', 'NFATC2')

sample_acts_mat <- tf_acts %>%
  tidyr::pivot_wider(id_cols = 'condition',
                     names_from = 'source',
                     values_from = 'score') %>%
  tibble::column_to_rownames('condition') %>%
  as.matrix()

tfs <- tf_acts %>%
  dplyr::group_by(source) %>%
  dplyr::summarise(std = stats::sd(score)) %>%
  dplyr::arrange(-abs(std)) %>%
  filter(source %in% reg) %>%
  dplyr::pull(source)

sample_acts_mat <- sample_acts_mat[,tfs]
sample_acts_mat <- scale(sample_acts_mat)


mat <- t(sample_acts_mat)

group_vec <- ifelse(grepl('0Gy', colnames(mat)), '0Gy','s')

cellorder = c('0Gy','s')
mat_group_mean <- sapply(cellorder, function(g) {
  rowMeans(mat[, group_vec == g, drop = FALSE], na.rm = TRUE)
})
colnames(mat_group_mean) <- cellorder

col_fun <- colorRamp2(c(-1, 0, 1), c("#313695", "white", "#A50026"))


pdf('figures/Fig4d_KHUcellline_TF_viper_heatmap.pdf', width = 4, height = 5)

Heatmap(
  mat_group_mean,
  name = 'TF acticity',
  col = col_fun,
  rect_gp = gpar(col = "white", lwd = 1),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 14),
  column_names_gp = gpar(fontsize = 16),
  column_names_rot = 90
)

dev.off()



### Extended Data Fig 4b - TF expression in our celllines
tf.exp = read.csv('data/4d/KHUcellline_TFexp_data.csv')

tf.exp$Celltype <- factor(tf.exp$Celltype, levels = c("WT",'PL',"SEN"))

comparisons_list <- list(c("SEN","WT"))

ggplot(tf.exp, aes(x = Celltype, y = value, fill = Celltype)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.7, width = 0.6) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.6,
               position = position_dodge(width = 0.6)) +
  stat_compare_means(comparisons = comparisons_list,
                     method = "t.test", label = "p.format", hide.ns = FALSE) +
  facet_wrap(~ gene, nrow = 2, scales = "free_y") +
  scale_fill_manual(values = c("WT" = "grey",  "SEN" = "tomato")) +
  labs(title = "Expression in KHU cell line",
       x = "Celltype", y = "Gene expression", fill = "Celltype") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 15))

ggsave('figures/FigSb_KHUcellline_TF_exp.pdf', width=6, height = 6, dpi=300)



