
setwd('/data_hdd4/yuseong/cancer_aging/analysis/h.github/Fig2')

rm(list=ls())


library(qs)
library(Seurat) 
library(ggplot2)  
library(lsa)        
library(pheatmap) 


### load data
pat = qread('data/patient_all_processed.qs')
pat$CellType = factor(pat$CellType, levels = c('MuscleCell', 'T.NKcell','PlasmaCell', 'Bcell', 'Myeloid', 'Fibroblast','Endothelial', 'Epithelial'))


### Fig 2a - UMAP of all celltypes
a = DimPlot(pat, reduction= "umap", label.size=8,
            cols = c('Epithelial' = '#E31A1C',
                     "Endothelial" = '#1F78B4',
                     "Fibroblast" = '#FF7F00',
                     'Bcell' = '#FB9A99',
                     "T.NKcell" = '#FDBF6F',
                     "Myeloid" = '#A6CEE3',
                     "PlasmaCell" = '#B2DF8A',
                     "MuscleCell" = '#33A02C'
            ))

a1 = a + theme_void() + 
  theme(legend.position = "right") +
  expand_limits(x = -17, y = -17) + 
  coord_cartesian(clip = "off") 

a1 = a1 + theme(
  legend.position = "right",
  legend.text = element_text(size = 12, face = "bold"), 
  legend.title = element_text(size = 14, face = "bold")
) +
  guides(color = guide_legend(override.aes = list(size = 8)),
         fill = guide_legend(override.aes = list(size = 8))) +
  expand_limits(x = -17, y = -17) + 
  coord_cartesian(clip = "off") +
  annotate("text", x = 15, y = 4, label = paste0("Total cells : ", format(nrow(pat@meta.data), big.mark = ",")), angle = 0, size = 6, fontface = "bold")
a1

ggsave('figures/Fig2a_patient_sc_all_umap_celltype.pdf', a1, width = 12, height = 12, dpi=300)




### Extended Data Fig 2a - expression of cell type marker genes 
FeaturePlot(pat, features = c('DCN', 'FCER1G', 'BANK1','JCHAIN','CD2','MYL1','KRT17', 'VWF'), ncol=3)

ggsave('figures/ExtendedFig2a_patient_sc_all_markers_featureplt.pdf', width = 10, height = 10)




### Fig 2b - magnified view of sub2 clusters
a = DimPlot(pat)
sel = CellSelector(a)

pat@meta.data$barcode= rownames(pat@meta.data)

sub.so = subset(pat, subset = (barcode %in% sel) & (CellType %in% c('Epithelial','Endothelial')))
DimPlot(sub.so, group.by = 'seurat_clusters')

sub.so$celltype_d = ifelse(sub.so$seurat_clusters %in% c('20','21'), 'Epithelial sub1', 
                           ifelse(sub.so$seurat_clusters %in% c('34'), 'Epithelial sub2', 
                                  ifelse(sub.so$seurat_clusters %in% c('33'), 'Endothelial sub2', 
                                         'Endothelial sub1')))

DimPlot(sub.so, group.by = 'celltype_d', 
        cols = c('Epithelial sub2' = 'firebrick2', 'Endothelial sub2' = 'black',
                 'Epithelial sub1' = '#fa8072', 'Endothelial sub1' = '#1F78B4'))+ theme_void()

ggsave('figures/Fig2b_patient_sc_all_umap_smse_closeup.pdf', width = 8, height = 6, dpi=300)




### Fig 2c - gene expression in epithelial/endothelial clusters
markers = c('LAMC2','KRT18',"KRT8", "PTHLH", ## epi
            "RAMP2",'CD34','CDH5','ECSCR', ## endo
            'FN1','S100A4', 'ZEB2',## emt
            'CCL3', 'CSF2', 'CXCL2','IL1B', ## sasp
            'MME','PPID','MMP2', 'MMP3' ## UV
)


DotPlot(sub.so, features = markers, group.by = 'celltype_d') +
  RotatedAxis() +
  scale_color_gradient2(
    low = "lightgrey", 
    mid = "indianred2", 
    high = "firebrick", 
    midpoint = 0.6 
  )

ggsave('figures/Fig2c_patient_sc_all_closeup_dotplt_markers.pdf', width = 10, height = 5)




### Fig 2d - cosine similarity between epithelial/endothelial clusters
avg_exp <- AverageExpression(sub.so, group.by = "celltype_d", slot = "data")$RNA

hvg <- VariableFeatures(pat)
avg_exp_hvg <- avg_exp[hvg, ]

cos_sim_mat <- lsa::cosine(as.matrix(avg_exp_hvg))
cos_sim_mat[1:4, 1:4]

palette1 <- colorRampPalette(c("blue", "white"))(33)   
palette2 <- colorRampPalette(c("white", "indianred1"))(33)  
palette3 <- colorRampPalette(c("indianred1", "red"))(34)    


final_palette <- c(palette1, palette2, palette3)

breaks <- c(
  seq(0, 0.2, length.out = 34),          
  seq(0.2, 0.9, length.out = 34)[-1],    
  seq(0.9, 1.0, length.out = 35)[-1]      
)


pdf('figures/Fig2c_patient_sc_all_umap_smse_closeup_corr.pdf', width = 5, height = 4)

pheatmap(cos_sim_mat, 
         name = "cosine\nsimilarity",  
         color = final_palette, 
         breaks = breaks
)

dev.off()



