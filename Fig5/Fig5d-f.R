
setwd('/data_hdd4/yuseong/cancer_aging/analysis/h.github/')

rm(list = ls())

library(ggplot2)
library(SpatialExperiment)
library(SummarizedExperiment)
library(spacexr)
library(qs)
library(MuSiC)
library(Seurat)
library(SeuratData)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(patchwork)
library(dplyr)


### load data
C1 = qread('Fig5/data/5d-f/Pre_spatial.qs')
C2 = qread('Fig5/data/5d-f/Post_spatial.qs')

### run RCTD
ref = qread('Fig2/data/2a-2d/patient_all_processed.qs')

ref@meta.data$rctd_ct = ref$celltype_f
table(ref$rctd_ct)

ref$SAMPLE = sapply(strsplit(ref$Sample, split = '-'), function(x) {strsplit(paste0(x[1:4], collapse='.'), split = '_')[[1]][1]})
table(ref$SAMPLE)

# create ref
cnt = as.data.frame(table(ref$condition, ref$rctd_ct))

## pre
ct1 = cnt[cnt$Var1 == 'Pre' & cnt$Freq > 25, ]$Var2
ref1 = subset(ref, subset = SAMPLE %in% c( "HyPR.HN.01.1", "HyPR.HN.02.1", "HyPR.HN.03.1", "HyPR.HN.04.1") & rctd_ct %in% ct1)
counts1 <- GetAssayData(ref1, assay = "RNA", slot = "counts")

cell_types1 <- ref1$rctd_ct
names(cell_types1) <- colnames(ref1)
unique(cell_types1)
cell_types1 = factor(cell_types1)

reference1 <- Reference(counts1, cell_types1)

## post
ct2 = cnt[cnt$Var1 == 'Post' & cnt$Freq > 25, ]$Var2
ref2 = subset(ref, subset = SAMPLE %in% c( "HyPR.HN.01.3", "HyPR.HN.02.3","HyPR.HN.04.3") & rctd_ct %in% ct2)

counts2 <- GetAssayData(ref2, assay = "RNA", slot = "counts")

cell_types2 <- ref2$rctd_ct
cell_types2 <- dplyr::case_when(cell_types2 == 'T/NKcell' ~ 'T.NKcell',
                                TRUE ~ cell_types2)
names(cell_types2) <- colnames(ref2)
unique(cell_types2)
cell_types2 = factor(cell_types2)

reference2 <- Reference(counts2, cell_types2)


# Create SpatialExperiment.
coords1 <- GetTissueCoordinates(C1)
coords1 <- data.frame(
  x = as.numeric(coords1[, 1]),
  y = as.numeric(coords1[, 2]),
  row.names = rownames(coords1)
)
counts1 <- GetAssayData(C1, assay = 'SCT', slot = 'counts')

spatialRNA1 <- SpatialRNA(coords1, counts1)


coords2 <- GetTissueCoordinates(C2)
coords2 <- data.frame(
  x = as.numeric(coords2[, 1]),
  y = as.numeric(coords2[, 2]),
  row.names = rownames(coords2)
)

counts2 <- GetAssayData(C2, assay = 'SCT', slot = 'counts')

spatialRNA2 <- SpatialRNA(coords2, counts2)


# preprocess data
rctd_data1 <- create.RCTD(spatialRNA1, reference1)
rctd_data2 <- create.RCTD(spatialRNA2, reference2)


# run RCTD
results_spe1 <- run.RCTD(rctd_data1,  doublet_mode = "doublet")
results_spe2 <- run.RCTD(rctd_data2,  doublet_mode = "doublet")

# add metadata
weight1 = data.frame(spacexr::normalize_weights(results_spe1@results$weights))
weight1 <- weight1[colnames(C1), ]
write.csv(weight1, 'Fig5/data/5d-f/Pre_RCTD_prop_normalized.csv', quote = F)

C1 <- AddMetaData(C1, metadata = weight1)
C1 <- AddMetaData(C1, metadata = results_spe1@results$results_df)


weight2 = data.frame(spacexr::normalize_weights(results_spe2@results$weights))
weight2 <- weight2[colnames(C2), ]
write.csv(weight2, 'Fig5/data/5d-f/Post_RCTD_prop_normalized.csv', quote = F)

C2 <- AddMetaData(C2, metadata = weight2)
C2 <- AddMetaData(C2, metadata = results_spe2@results$results_df)


qsave(C1, 'Fig5/data/5d-f/Pre_spatial_RCTD+prop.qs', preset='fast')
qsave(C2, 'Fig5/data/5d-f/Post_spatial_RCTD+prop.qs', preset='fast')



#### Fig 5d - RCTD results in post sample
SpatialDimPlot(C2, group.by = 'first_type', pt.size=4) &
  scale_fill_manual(
    values = list(
      'Bcell' = '#FB9A99',
      "T.NKcell" = '#FDBF6F',
      "Myeloid" = '#A6CEE3',
      "Fibroblast" = '#33A02C',
      "PlasmaCell" = '#B2DF8A',
      "Endothelial" = '#1F78B4',
      "MuscleCell" = '#FF7F00',
      'Sen_Endothelial' = '#000000',
      'Sen_Malignant' = '#FF3030',
      'Non_Malignant' = '#FFE4C4',
      'Malignant' = '#6A5ACD'))

ggsave('Fig5/figures/Fig5d_Post_RCTD_celltypes_first_type.pdf', width=12, height=8, dpi=300)




#### Fig 5e - RCTD results of senescent cells in post sample
SpatialFeaturePlot(C2, features = "Sen_Malignant", pt.size.factor = 5.5) +
  scale_fill_gradientn(colors = c("gray95",  "orange", "red", "darkred"))

ggsave('Fig5/figures/Fig5e_Post_RCTD_SenMal.pdf', width=8, height=6, dpi=300)


SpatialFeaturePlot(C2, features = "Sen_Endothelial", pt.size.factor = 5.5) +
  scale_fill_gradientn(colors = c("gray95", "gray80", "gray60", "gray40", "black"))

ggsave('Fig5/figures/Fig5e_Post_RCTD_SenEndo.pdf', width=8, height=6, dpi=300)




#### Fig 5f - correlation between proportion of celltypes per spot in post sample
post_weight = read.csv('Fig5/data/5d-f/Post_RCTD_prop_normalized.csv', row.names = 1)

cor_mat = cor(post_weight, method = "pearson")

ggcorrplot(cor_mat, hc.order = F, outline.color = "white", tl.srt = 60,
           tl.cex = 18, lab_size = 7, colors = c("#91A28C","white","#8F2C37")) +
  theme(plot.margin = ggplot2::margin(0.1, 0.1, 0.1, 0.1, "cm"), panel.background = element_blank(),
        plot.background = element_blank(), panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
        axis.text = element_text(size = 12),  axis.ticks =element_blank(), axis.title =element_blank(),
        legend.title=element_text(size = 16,face="bold"), legend.text=element_text(size = 16),
        legend.key = element_rect(colour = "transparent", fill = "white"), legend.key.size = unit(0.45, 'cm')) +
  coord_fixed() + ggtitle("Correlation") + theme(plot.title = element_text(size=22,face="bold"))

ggsave('Fig5/figures/Fig5f_Post_cell_colocalization_Normalized.pdf', width=8, height =6, dpi=300)




### Extended Data Fig 5a - RCTD results per celltype in pre sample
for (ct in c("Bcell","Endothelial","Fibroblast",
             "Malignant","MuscleCell","Myeloid",
             "Non_Malignant","PlasmaCell","T.NKcell")) {
  SpatialFeaturePlot(C1, features = ct, pt.size.factor = 5.5)
  ggsave(paste0('Fig5/figures/FigS5a_Pre_RCTD_',ct,'.pdf'), width=8, height=6, dpi=300)
}




### Extended Data Fig 5b - RCTD results in pre sample
SpatialDimPlot(C1, group.by = 'first_type', pt.size=4) &
  scale_fill_manual(
    values = list(
      'Bcell' = '#FB9A99',
      "T.NKcell" = '#FDBF6F',
      "Myeloid" = '#A6CEE3',
      "Fibroblast" = '#33A02C',
      "PlasmaCell" = '#B2DF8A',
      "Endothelial" = '#1F78B4',
      "MuscleCell" = '#FF7F00',
      'Sen_Endothelial' = '#000000',
      'Sen_Malignant' = '#FF3030',
      'Non_Malignant' = '#FFE4C4',
      'Malignant' = '#6A5ACD'))

ggsave('Fig5/figures/FigS5b_RCTD_celltypes_first_type.pdf', width=12, height=8, dpi=300)




### Extended Data Fig 5c - RCTD results per celltype in post sample
for (ct in c("Bcell","Endothelial","Fibroblast",
             "Malignant","MuscleCell","Myeloid",
             "Non_Malignant","PlasmaCell","T.NKcell")) {
  SpatialFeaturePlot(C2, features = ct, pt.size.factor = 5.5)
  ggsave(paste0('Fig5/figures/FigS5c_Post_RCTD_',ct,'.pdf'), width=8, height=6, dpi=300)
}

