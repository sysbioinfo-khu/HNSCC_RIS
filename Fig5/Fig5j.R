setwd('/data_hdd4/yuseong/cancer_aging/analysis/h.github/Fig5/')

rm(list=ls())

library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(NICHES)
library(dplyr)
library(ggplot2)
library(qs)
library(tidyr)


### run NICHES
seurat_object = qread('data/5j/ST_decompositon_5_seuratobject.qs')

mal.near = read.csv('data/5j/spotiphy_malignant_near_info.csv', row.names = 1)
mal.near$near3 = ifelse(mal.near$near3 == 'exclude', 'exclude_mal', mal.near$near3)
end.near = read.csv('data/5j/spotiphy_endothelial_near_info.csv', row.names = 1)
end.near$near3 = ifelse(end.near$near3 == 'exclude', 'exclude_end', end.near$near3)
head(end.near)


seurat_object@meta.data$cell_num = rownames(seurat_object@meta.data)

seurat_object$cell_type <- as.character(seurat_object$cell_type)

seurat_object@meta.data <- seurat_object@meta.data %>%
  mutate(cell_type_n = case_when(
    cell_num %in% rownames(mal.near) ~ mal.near[cell_num, "near3"],
    cell_num %in% rownames(end.near) ~ end.near[cell_num, "near3"],
    TRUE ~ cell_type  
  ))
table(seurat_object$cell_type_n)
table(mal.near$near3)                                                
table(end.near$near3)                                                

Idents(seurat_object) <- seurat_object$cell_type_n
sub <- subset(seurat_object,idents = c('Sen_Malignant','Sen_Endothelial','near_mal','near_endo','far_mal','far_endo'))
DimPlot(sub,reduction = 'umap',label = T,repel = F,label.size = 6)


lr_network = read.csv('data/5j/merged_lr_network_for_nichenetr.csv')

lr.dat = as.data.frame(lr_network)
colnames(lr.dat) = c('ligand','receptor')
rownames(lr.dat) = paste0(lr.dat$ligand, '_', lr.dat$receptor)
head(lr.dat)


sub$celltype = ifelse(sub$cell_type_n %in% c('Sen_Malignant','Sen_Endothelial'), 'SEN', sub$cell_type_n)

imputed <- SeuratWrappers::RunALRA(sub)
scc.imputed <- RunNICHES(imputed,
                         assay = 'alra',
                         species = 'human',
                         LR.database = 'custom',
                         cell_types = 'celltype',
                         custom_LR_database = lr.dat,
                         CellToCell = T)
demo.2 <- scc.imputed$CellToCell
demo.2 <- ScaleData(demo.2)
demo.2 <- RunPCA(demo.2,features = rownames(demo.2))
ElbowPlot(demo.2,ndims=50)

PCHeatmap(demo.2,dims = 1:10,balanced = T,cells = 100)

demo.2 <- RunUMAP(demo.2,dims = 1:10)
DimPlot(demo.2,reduction = 'umap',group.by = 'VectorType',label = F)

qsave(demo.2, 'data/5j/runNiches_demo_imputed_merged.qs', preset='fast')



### visualization
demo.2 = qread('data/5j/runNiches_demo_imputed_merged.qs')

# Find markers
mark <- FindAllMarkers(demo.2,min.pct = 0.25,only.pos = T,test.use = "roc")
GOI_niche <- mark %>% group_by(cluster) %>% top_n(5,myAUC)

demo.subset <- subset(demo.2, idents = c("SEN—near_endo", "SEN—far_endo", 
                                         "SEN—near_mal", "SEN—far_mal"))

demo.subset <- ScaleData(demo.subset)
DoHeatmap(demo.subset, features = unique(GOI_niche$gene)) +
  scale_fill_gradientn(colors = c("grey", "white", "blue"))

## interactions via NichenetR
sel.m.p = c('CTHRC1—FZD5','CTHRC1—FZD3','CTHRC1—ROR2',
            'FN1—C5AR1', 'FN1—ITGA4',
            'FN1—ITGA5','FN1—ITGB7','FN1—ITGB8','FN1—NT5E',
            'FN1—PLAUR','FN1—TNFRSF11B', 'FN1—ITGB1','FN1—ITGA2',
            'MMP2—EGFR', 'MMP2—DDR1','MMP2—FGFR1','MMP2—PECAM1','MMP2—F2R', 
            'TIMP1—FGFR2','TIMP1—MMP1')

sel.e.p = c('FN1—ITGAV','FN1—SDC2','FN1—CD79A','FN1—COL13A1',
            'FN1—IL17RC','FN1—ITGA6','FN1—ITGA8','FN1—ITGA9',
            'FN1—NT5E','FN1—PLAUR','FN1—ROBO4','FN1—TNFRSF11B','FN1—SDC1',
            'MMP1—F2R', 'MMP1—CD44',
            'MMP2—F2R','MMP2—ITGAV','MMP2—SDC2',
            'TIMP1—MMP1','TIMP1—MMP9')

sel.pairs.in = intersect(sel.m.p, sel.e.p)



plot_data <- DotPlot(demo.subset, features = sel.pairs.in)$data

### endo
plot_data.endo <- plot_data %>% filter(grepl('endo', id) & features.plot %in% sel.pairs.in) %>% 
  separate(id, into = c("Cell_Group", "Location"), sep = "—", remove = FALSE) %>%
  mutate(`SendRe` = case_when(
    grepl('mal', Location)  ~ paste0(Cell_Group, ' - Mal'),
    grepl('endo', Location) ~ paste0(Cell_Group, ' - Endo'),
    TRUE                    ~ Cell_Group
  ))
unique(plot_data.endo$SendRe)

plot_data.endo$Location = factor(plot_data.endo$Location, levels = c('near_endo','far_endo'))
p <- ggplot(plot_data.endo, aes(x = features.plot, y = Location, size = pct.exp, color = avg.exp.scaled)) +
  geom_point() +
  scale_color_gradient2(
    low = 'blue', 
    mid = 'white', 
    high = 'red', 
    midpoint = 0,
    limits = c(-1, 1),
    oob = scales::squish
  ) +
  scale_size_continuous(
    range = c(1, 6),
    limits = c(0, 100),                      
    breaks = c(25, 50, 75, 100),          
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold")) +
  labs(title = "NEAR ENDO Interaction",
       x = "Ligand-Receptor Pairs",
       y = "Location",
       color = "Avg Expression",
       size = "Percent Expressed") +
  coord_flip() +
  facet_wrap(~SendRe, scales = "free_y", ncol = 4)

p


### mal
plot_data.mal <- plot_data %>% filter(grepl('mal', id) & features.plot %in% sel.pairs.in) %>% 
  separate(id, into = c("Cell_Group", "Location"), sep = "—", remove = FALSE) %>%
  mutate(`SendRe` = case_when(
    grepl('mal', Location)  ~ paste0(Cell_Group, ' - Mal'),
    grepl('endo', Location) ~ paste0(Cell_Group, ' - Endo'),
    TRUE                    ~ Cell_Group
  ))
unique(plot_data.mal$SendRe)


plot_data.mal$Location = factor(plot_data.mal$Location, levels = c('near_mal','far_mal'))

p2 <- ggplot(plot_data.mal, aes(x = features.plot, y = Location, size = pct.exp, color = avg.exp.scaled)) +
  geom_point() +
  scale_color_gradient2(
    low = 'blue', 
    mid = 'white', 
    high = 'red', 
    midpoint = 0,
    limits = c(-1, 1),
    oob = scales::squish
  ) +
  scale_size_continuous(
    range = c(1, 6),
    limits = c(0, 100),                      
    breaks = c(25, 50, 75, 100),          
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold")) +
  labs(title = "NEAR MAL Interactions",
       x = "Ligand-Receptor Pairs",
       y = "Location",
       color = "Avg Expression",
       size = "Percent Expressed") +
  coord_flip() +
  facet_wrap(~SendRe, scales = "free_y", ncol = 4)
p2

p2+p

ggsave('figures/Fig5j_niches_plot_merged.pdf', p2+p, width=10, height=5, dpi=300)




