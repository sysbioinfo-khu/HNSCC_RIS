library(dplyr)
library(tidyr)
library(ggplot2)
library(limma)
library(Biobase)
library(umap)
library(ggrepel)
library(sigminer)

rm(list=ls())

setwd('/data_hdd4/yuseong/cancer_aging/analysis/h.github/Fig3/')


### extract HyHR.PN.02 clone D specific CNVs
ol.cnv.df2 = read.csv('data/3c/clone_analysis/ALL_sample_clones_cnv_results_+DNAmicroarray.csv')
col2 = colnames(ol.cnv.df2)[grepl('HyHR.PN.02', colnames(ol.cnv.df2))]
col2 = c("HyHR.PN.02_Clone.A",
         "HyHR.PN.02_Clone.B",
         "HyHR.PN.02_Clone.C",
         "HyHR.PN.02_Clone.D",
         "HyHR.PN.02_Clone.E",
         "HyHR.PN.02_Clone.F",
         "HyHR.PN.02_Clone.A.HyHR.PN.02_Clone.B.HyHR.PN.02_Clone.C.HyHR.PN.02_Clone.D.HyHR.PN.02_Clone.E.HyHR.PN.02_Clone.F",
         "HyHR.PN.02_Clone.B.HyHR.PN.02_Clone.C.HyHR.PN.02_Clone.D.HyHR.PN.02_Clone.E.HyHR.PN.02_Clone.F",
         "HyHR.PN.02_Clone.B.HyHR.PN.02_Clone.C",
         "HyHR.PN.02_Clone.D.HyHR.PN.02_Clone.E.HyHR.PN.02_Clone.F",
         "HyHR.PN.02_Clone.E.HyHR.PN.02_Clone.F")

ol.cnv.df2 = ol.cnv.df2[colnames(ol.cnv.df2) %in% col2]
cnv2 = c()
for (col in colnames(ol.cnv.df2)) {
  cnv2 = c(cnv2, unique(ol.cnv.df2[[col]]))
}
cnv2 = na.omit(cnv2)

cld.cnv = na.omit(ol.cnv.df2["HyHR.PN.02_Clone.D"])
cld.cnv = unique(cld.cnv$HyHR.PN.02_Clone.D)

### load DNAmicroarray DE results
tT_b = read.csv('data/S3c/DNAmicroarray_limma_results.csv')
head(tT_b)


### plot DE results and HyHR.PN.02 clone D specific CNVs
tT_b <- tT_b %>%
  filter(!grepl('X|Y', LIFTOVER_CHROMOSOME)) %>%
  mutate(
    band_cnv = paste0(LIFTOVER_CYTOBAND,diffexpressed),
    Annotation = ifelse(band_cnv %in% cnv2, band_cnv, NA)
  )

tT_b$diffexpressed = ifelse(((tT_b$diffexpressed == 'gain') & (tT_b$band_cnv %in% cld.cnv)), 'cloneD_gain',
                            ifelse(((tT_b$diffexpressed == 'loss') & (tT_b$band_cnv %in% cld.cnv)), 'cloneD_loss',
                                   tT_b$diffexpressed))

table(tT_b$diffexpressed)

tT_f = tT_b[tT_b$diffexpressed != 'NO',]

ggplot(data = tT_b, aes(x = logFC, y = -log10(P.Value), color = diffexpressed, label = LIFTOVER_CYTOBAND)) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey") +
  xlim(-2,2)+
  scale_color_manual(values = c("cloneD_loss" = "red", "cloneD_gain" = "blue", 
                                'loss' = 'lightpink2', 'gain' = 'lightskyblue', "NS" = "grey")) +
  labs(
    x = expression("log"[2]*"(fold change)"),
    y = expression("-log"[10]*"(p-value)"),
    color = "CNV Status"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave('figures/DE/DNAmicroarray_volcanoplot.pdf', width=8, height = 6, dpi=300)

