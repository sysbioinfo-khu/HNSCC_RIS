
setwd('/data_hdd4/yuseong/cancer_aging/analysis/h.github/Fig6/')

rm(list=ls())

library(dplyr)
library(qs)
library(Seurat)
library(SeuratObject)
library(patchwork)
library(escape)
library(BiocParallel)
library(ggpubr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(readxl)
library(ComplexHeatmap)
library(circlize)
library(pals)
library(viridis)

fib.a = qread('data/6a-c/patient_seurat_fibroblast_subclusters.qs')


#######
jaccard_index <- function(set1, set2) {
  inter <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  if (union == 0) return(0)
  inter / union
}


cluster_markers = readRDS('data/6a-c/cluster_markers_200.rds')
caf_markers = readRDS('data/6a-c//caf_markers_200.rds')

mat <- matrix(0, nrow = length(cluster_markers), ncol = length(caf_markers))
rownames(mat) <- names(cluster_markers)
colnames(mat) <- names(caf_markers)

for (i in seq_along(cluster_markers)) {
  for (j in seq_along(caf_markers)) {
    mat[i, j] <- jaccard_index(cluster_markers[[i]], caf_markers[[j]])
  }
}

rng <- range(mat, na.rm = TRUE)
mat01 <- if (diff(rng) == 0) {
  matrix(0, nrow(mat), ncol(mat), dimnames = dimnames(mat))
} else {
  (mat - rng[1]) / diff(rng)
}

mat01 = mat01[c('F0','F1','F2','F3','F4'), c('Interferon.response.CAF',
                                             'Myofibroblast',
                                             'Progenitor.like.Fib',
                                             'Inflammatory.CAF',
                                             'Antigen.presenting.CAF',
                                             'Vascular.CAF'
)]

pal <- colorRampPalette(rev(brewer.pal(11, "Spectral")))(200)


pdf('figures/Fig6b_patient_sc_fib_clusters_jaccard_index_spectral_200.pdf', width=8, height = 4)

h = Heatmap(mat01,
            col = pal,
            name = "Scaled\nJaccard\nIndex",
            row_names_side = "left",
            column_names_side = "bottom",
            cluster_rows = F,
            cluster_columns = F,
            # col = colorRamp2(c(0, 0.4, 0.6,1), c("skyblue", "white", "lightsalmon1", 'sienna1')),
            column_names_rot = 45,
)

draw(h, padding = unit(c(3, 15, 3, 3), "mm"))

dev.off()


###################
############################
## FIBROBLAST REF MARKERS ##
############################
##----------------------------------------
## cross-tissue human fibroblast ~~
# Gao, Yang et al. 
# “Cross-tissue human fibroblast atlas reveals myofibroblast subtypes with distinct roles in immune modulation.”
# Cancer cell vol. 42,10 (2024): 1764-1783.e10. doi:10.1016/j.ccell.2024.08.020
# S2 Progenitor.like.Fibroblast
# S4 Inflammatory.Fibroblast
## ---------------------------------------

#### 
caf_markers['Progenitor.like.Fib'] = caf_markers['c03'] # COL15A1+Progenitor.like.fib
caf_markers['Myofibroblast'] = caf_markers['c04'] # LRRC15+Myofib
caf_markers['Inflammatory.CAF'] = caf_markers['c07'] # IL6+Inflammatory.fib


##---------------------------------------
## breast cancer
# Cords, L., Tietscher, S., Anzeneder, T. et al.
# Cancer-associated fibroblast classification in single-cell and spatial proteomics data.
# Nat Commun 14, 4294 (2023). https://doi.org/10.1038/s41467-023-39762-1
# S3 apCAF
# S4 tCAF
# S6 vCAF
# S7 dCAF
##---------------------------------------

### 
caf_markers['Antigen.presenting.CAF'] = caf_markers['apCAF']
caf_markers['Vascular.CAF'] = caf_markers['vCAF']
caf_markers['Interferon.response.CAF'] = caf_markers['ifnCAF']

