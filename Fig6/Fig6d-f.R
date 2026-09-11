
setwd('/data_hdd4/yuseong/cancer_aging/analysis/h.github/')

rm(list=ls())


library(qs)    
library(Seurat)       
library(dplyr)     
library(tidyr)
library(tibble)
library(SPATA2)       
library(ggplot2)     



## load data
fib.so = qread('Fig6/data/6d-f/spotiphy_fibroblast_processed.qs')
fib.sc = qread('Fig6/data/6a-c/patient_seurat_fibroblast_subclusters.qs')


## calculate proportion of fibroblast subtypes in spot
anchors <- FindTransferAnchors(reference = fib.sc, query = fib.so, dims = 1:24)

fib.pred <- TransferData(anchorset = anchors, refdata = fib.sc$Fibroblast_sub, 
                         dims = 1:24)
table(fib.pred$predicted.id)

fib.so <- AddMetaData(fib.so, metadata = fib.pred)

fib.cols = list(
  'F0' = '#EEDD82',
  'F1' ='#96CDCD',
  'F2' = '#EE799F',
  'F3' = '#66CDAA',
  'F4' = '#83B366'
)


DimPlot(fib.so, group.by = 'predicted.id', cols = fib.cols) + theme_void() + 
  theme(legend.position = "right") +
  coord_cartesian(clip = "off") 

table(fib.so$predicted.id)

fib.spot = fib.so@meta.data %>% group_by(predicted.id, spot_name) %>% summarise(fib.cnt = n()) %>%
  group_by(spot_name) %>%
  slice_max(fib.cnt, n = 1, with_ties = FALSE) %>% 
  ungroup()


C2 = qread('Fig5/data/5d-f/Post_spatial_RCTD+prop.qs')
C2@meta.data$barcode = rownames(C2@meta.data)

C2@meta.data$fib.anno = ifelse(C2$barcode %in% unique(fib.spot[fib.spot$predicted.id == 'F0',]$spot_name), 'F0',
                               ifelse(C2$barcode %in% unique(fib.spot[fib.spot$predicted.id == 'F1',]$spot_name), 'F1',
                                      ifelse(C2$barcode %in% unique(fib.spot[fib.spot$predicted.id == 'F2',]$spot_name), 'F2',
                                             ifelse(C2$barcode %in% unique(fib.spot[fib.spot$predicted.id == 'F3',]$spot_name), 'F3',
                                                    ifelse(C2$barcode %in% unique(fib.spot[fib.spot$predicted.id == 'F4',]$spot_name), 'F4',
                                                           'none')))))

table(C2@meta.data$fib.anno)

C2 <- UpdateSeuratObject(C2)

SpatialDimPlot(C2, group.by = 'fib.anno', pt.size.factor=4, cols = fib.cols)

fib.prop = fib.so@meta.data %>% group_by(spot_name) %>% 
  summarise(cnt = n(),
            F0 = sum(prediction.score.F0) / cnt,
            F1 = sum(prediction.score.F1) / cnt,
            F2 = sum(prediction.score.F2) / cnt,
            F3 = sum(prediction.score.F3) / cnt,
            F4 = sum(prediction.score.F4) / cnt,
  ) %>% 
  mutate(total_sum = F0 + F1 + F2 + F3 + F4,
         F0 = F0 / total_sum,
         F1 = F1 / total_sum,
         F2 = F2 / total_sum,
         F3 = F3 / total_sum,
         F4 = F4 / total_sum
  ) %>% dplyr::select(-c('total_sum')) %>% as.data.frame()
rownames(fib.prop) = fib.prop$spot_name
fib.prop$spot_name = NULL

head(fib.prop)

C2 = AddMetaData(C2, fib.prop)
colnames(C2@meta.data)

C2$fib.F0 = C2$F0 * C2$Fibroblast
C2$fib.F1 = C2$F1 * C2$Fibroblast
C2$fib.F2 = C2$F2 * C2$Fibroblast
C2$fib.F3 = C2$F3 * C2$Fibroblast
C2$fib.F4 = C2$F4 * C2$Fibroblast

C2@meta.data$celltype = C2@meta.data$first_type
table(C2@meta.data$celltype)


#### Fig 6d - fibroblast subtypes in 10x visium post sample
SpatialDimPlot(C2, group.by = 'fib.anno', pt.size.factor = 5.5, cols = c('F0' = '#EEDD82',
                                                                         'F1' ='#96CDCD',
                                                                         'F2' = '#EE799F',
                                                                         'F3' = '#66CDAA',
                                                                         'F4' = '#83B366'
)) +
  guides(color = guide_legend(override.aes = list(size = 8)),
         fill = guide_legend(override.aes = list(size = 8)))

ggsave('Fig6/figures/Fig6d_Post_FIBsubtype_dimplot.pdf', width=6, height = 8, dpi=300)

qsave(C2, 'Fig6/data/6d-f/Post_spatial_RCTD+prop_+spotiphy_fib.anno.qs', preset='fast')


## distance of fibroblast subtypes from sen niche
C2 = qread('Fig6/data/6d-f/Post_spatial_RCTD+prop_+spotiphy_fib.anno.qs')

C2_spata = asSPATA2(object = C2, sample_name = 'C2', assay_modality = "gene", img_scale_fct = 'lowres', assay_name = 'Spatial', platform = "VisiumSmall")
rownames(C2_spata@meta_obs) = NULL

C2_spata <- rotateAll(C2_spata, angle = 90, clockwise = TRUE)
C2_spata <- flipCoordinates(C2_spata, axis = "vertical")

plotSurface(C2_spata, color_by = "celltype", pt_size = 1.2)
plotSurface(C2_spata, color_by = "fib.anno", pt_clrp = "npg")

table(C2_spata@meta_obs$fib.anno)



C2_spata = createGroupAnnotations(
  object = C2_spata, 
  grouping = "celltype", 
  group = c("Sen_Malignant", "Sen_Endothelial"),
  id = "Sen", 
  overwrite = TRUE
)

getSpatAnnIds(C2_spata)

# upper region
senmal_outline <- 
  ggpLayerSpatAnnOutline(C2_spata, ids = "Sen_1")
plotImage(C2_spata) + 
  senmal_outline

plotSurface(C2_spata, color_by = "Sen_Malignant", pt_clrsp = "plasma") + 
  senmal_outline
plotSurface(C2_spata, color_by = "Sen_Endothelial", pt_clrsp = "plasma") + 
  senmal_outline


C2_spata = identifyTissueOutline(C2_spata, method = "obs")

coords_df <- 
  getCoordsDfSA(C2_spata, ids = c('Sen_1'), unit = "mm", binwidth = "200um", core0 = T) 


C2_spata <- 
  addFeatures(C2_spata, feature_df = coords_df[,c("barcodes", "dist", "bins_dist")], overwrite = T)

sas_screning <- 
  ggpLayerScreeningDirectionSAS(C2_spata, ids = c('Sen_1'), line_size = 0.75)


plotImage(C2_spata) + 
  senmal_outline + 
  sas_screning

plotSurface(C2_spata, color_by = "dist", pt_clrsp = "plasma") + 
  senmal_outline + 
  sas_screning + 
  labs(color = "Dist [mm]")


#### Fig 6e left - sen malignant proportion in sen niche
plotSurface(C2_spata, color_by = "Sen_Malignant", display_image = F, outline = T, pt_clrsp = "BuPu") + 
  senmal_outline + 
  sas_screning

ggsave('Fig6/figures/Fig6e_Post_sen_endo_sen_mal_distance_plot.pdf', width=6, height = 8, dpi=300)


#### Fig 6e left - sen malignant proportion in sen niche
plotSurface(C2_spata, color_by = "Sen_Endothelial", display_image = F, outline = T, pt_clrsp = "BuPu") + 
  senmal_outline + 
  sas_screning

ggsave('Fig6/figures/Fig6e_Post_sen_endo_distance_plot.pdf', width=6, height = 8, dpi=300)




#### Fig 6f - fibroblast proportion according to the distance from sen niche
plotSasLineplot(C2_spata, ids = c('Sen_1'), variables = "fib.F0", line_color = "lightgoldenrod") + 
  labs(x = "Distance to Sen Borders [mm]") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave('Fig6/figures/Fig6f_Post_sen_mal_distance_plot_S1_nopanel.pdf', width=6, height = 4, dpi=300)

plotSasLineplot(C2_spata, ids = c('Sen_1'), variables = "fib.F1", line_color = "paleturquoise3") + 
  labs(x = "Distance to Sen Borders [mm]")+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave('Fig6/figures/Fig6f_Post_sen_mal_distance_plot_S2_nopanel.pdf', width=6, height = 4, dpi=300)

plotSasLineplot(C2_spata, ids = c('Sen_1'), variables = "fib.F2", line_color = "palevioletred2") + 
  labs(x = "Distance to Sen Borders [mm]")+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave('Fig6/figures/Fig6f_Post_sen_mal_distance_plot_S3_nopanel.pdf', width=6, height = 4, dpi=300)

plotSasLineplot(C2_spata, ids = c('Sen_1'), variables = "fib.F3", line_color = "mediumaquamarine") + 
  labs(x = "Distance to Sen Borders [mm]")+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave('Fig6/figures/Fig6f_Post_sen_mal_distance_plot_S4_nopanel.pdf', width=6, height = 4, dpi=300)

plotSasLineplot(C2_spata, ids = c('Sen_1'), variables = "fib.F4", line_color = "#83B366") + 
  labs(x = "Distance to Sen Borders [mm]")+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave('Fig6/figures/Fig6f_Post_sen_mal_distance_plot_S5_nopanel.pdf', width=6, height = 4, dpi=300)




#### Fig 6g - iCAF marker gene expression according to the distance from sen niche
plotSasLineplot(C2_spata, ids = c('Sen_1'), variables = "IL6", line_color = "palevioletred2") + 
  labs(x = "Distance to Sen Borders [mm]")+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave('Fig6/figures/Fig6g_Post_sen_mal_distance_plot_IL6_nopanel.pdf', width=6, height = 4, dpi=300)

plotSasLineplot(C2_spata, ids = c('Sen_1'), variables = "CXCL1", line_color = "palevioletred2") + 
  labs(x = "Distance to Sen Borders [mm]")+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave('Fig6/figures/Fig6g_Post_sen_mal_distance_plot_CXCL1_nopanel.pdf', width=6, height = 4, dpi=300)



# down region
senmal_outline <- 
  ggpLayerSpatAnnOutline(C2_spata, ids = "Sen_2")
plotImage(C2_spata) + 
  senmal_outline

plotSurface(C2_spata, color_by = "Sen_Malignant", pt_clrsp = "plasma") + 
  senmal_outline
plotSurface(C2_spata, color_by = "Sen_Endothelial", pt_clrsp = "plasma") + 
  senmal_outline


C2_spata = identifyTissueOutline(C2_spata, method = "obs")

coords_df <- 
  getCoordsDfSA(C2_spata, ids = c('Sen_2'), unit = "mm", binwidth = "200um", core0 = T) 

C2_spata <- 
  addFeatures(C2_spata, feature_df = coords_df[,c("barcodes", "dist", "bins_dist")], overwrite = T)

sas_screning <- 
  ggpLayerScreeningDirectionSAS(C2_spata, ids = c('Sen_2'), line_size = 0.75)


plotImage(C2_spata) + 
  senmal_outline + 
  sas_screning

plotSurface(C2_spata, color_by = "dist", pt_clrsp = "plasma") + 
  senmal_outline + 
  sas_screning + 
  labs(color = "Dist [mm]")


#### Extended Data Fig 6f - sen malignant proportion in sen niche2
plotSurface(C2_spata, color_by = "Sen_Malignant", display_image = F, outline = T, pt_clrsp = "BuPu") + 
  senmal_outline + 
  sas_screning

ggsave('Fig6/figures/FigS6f_Post_sen_endo_sen_mal_distance_plot_sen2.pdf', width=6, height = 8, dpi=300)


#### Extended Data Fig 6f - sen endothelial proportion in sen niche2
plotSurface(C2_spata, color_by = "Sen_Endothelial", display_image = F, outline = T, pt_clrsp = "BuPu") + 
  senmal_outline + 
  sas_screning

ggsave('Fig6/figures/FigS6f_Post_sen_endo_distance_plot_sen2.pdf', width=6, height = 8, dpi=300)




#### Extended Data Fig 6g - fibroblast proportion difference according to the distance from sen niche2
plotSasLineplot(C2_spata, ids = c('Sen_2'), variables = "fib.F2", line_color = "palevioletred2") + 
  labs(x = "Distance to Sen Borders [mm]")+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave('Fig6/figures/FigS6g_Post_sen_mal_distance_plot_S3_nopanel_sen2.pdf', width=6, height = 4, dpi=300)


#### Extended Data Fig 6h - iCAF marker gene expression according to the distance from sen niche2
plotSasLineplot(C2_spata, ids = c('Sen_2'), variables = "IL6", line_color = "palevioletred2") + 
  labs(x = "Distance to Sen Borders [mm]")+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave('Fig6/figures/FigS6h_Post_sen_mal_distance_plot_IL6_nopanel_sen2.pdf', width=6, height = 4, dpi=300)

plotSasLineplot(C2_spata, ids = c('Sen_2'), variables = "CXCL1", line_color = "palevioletred2") + 
  labs(x = "Distance to Sen Borders [mm]")+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave('Fig6/figures/FigS6h_Post_sen_mal_distance_plot_CXCL1_nopanel_sen2.pdf', width=6, height = 4, dpi=300)
