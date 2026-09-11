setwd('/data_hdd4/yuseong/cancer_aging/analysis/h.github')

rm(list=ls())


library(qs)
library(Seurat)
library(ggpubr)


mal.so = qread('Fig2/data/2e_g_i/patient_mal_seurat_+GSVA+ME_ssgsea_bicor.qs')

### GENE LIST
## Noh, Joo Kyung et al.
## “SOD2- and NRF2-associated Gene Signature to Predict Radioresistance 
## in Head and Neck Cancer.” 
## Cancer genomics & proteomics vol. 18,5 (2021): 675-684. doi:10.21873/cgp.20289

rr_scores = read.csv('Fig3/data/3a/patient_mal_ssGSEA_scores_radioresistance_norm.csv', row.names=1)

mal.so = AddMetaData(mal.so, rr_scores)

mal.so$patient = sapply(strsplit(mal.so$Sample, split = '-'), function(x) {paste0(x[1:3], collapse='-')})
table(mal.so$patient)

comp_p = list(c("HyPR-HN-01", "HyPR-HN-02"), 
              c("HyPR-HN-03", "HyPR-HN-02"),
              c("HyPR-HN-04", "HyPR-HN-02"))

VlnPlot(subset(mal.so, subset = condition == 'Pre'), features = 'Radioresistance', group.by = 'patient', pt.size=0, cols = c(
  "HyPR-HN-01" = "#79CDCD", "HyPR-HN-02" = "#B22222", "HyPR-HN-03" = "#BEBEBE", "HyPR-HN-04" = "#8FBC8F"
)) &
  stat_summary(geom='crossbar', fun=mean, width=0.3) &
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) &
  stat_compare_means(comparisons = comp_p, method = "t.test", label = "p.signif")

ggsave('Fig3/figures/Fig3a_sc_mal_radioresistance_vlnplt_norm.pdf', width=6, height=4)


