setwd('/data_hdd4/yuseong/cancer_aging/analysis/h.github/')

rm(list=ls())

library(nichenetr)
library(RColorBrewer)
library(tidyverse)
library(Seurat) 
library(qs)
library(dplyr)
library(patchwork)
library(ggpubr)
library(cowplot)
library(scales)
library(harmony)
library(pheatmap)


#### load data
pat = qread('Fig2/data/2a-2d/patient_all_processed.qs')
DimPlot(pat, group.by='celltype_f', label=T)

pat[['cci_celltype']] = case_when(as.character(pat$celltype_f) == 'sen_malignant' ~ 'Sen_Malignant',
                                  as.character(pat$celltype_f) == 'sen_endothelial' ~ 'Sen_Endothelial',
                                  as.character(pat$celltype_f) == 'malignant' & pat$condition == 'Post' ~ 'Post_Malignant',
                                  as.character(pat$celltype_f) == 'malignant' & pat$condition == 'On' ~ 'On_Malignant',
                                  as.character(pat$celltype_f) == 'malignant' & pat$condition == 'Pre' ~ 'Pre_Malignant',
                                  as.character(pat$celltype_f) == 'Non_Malignant' ~ 'Non_Malignant',
                                  TRUE ~ as.character(pat$celltype_f))

table(pat$cci_celltype)

DimPlot(pat, group.by = 'cci_celltype')

pat = SetIdent(pat, value = "cci_celltype")

seuratObj <- alias_to_symbol_seurat(pat, "human")
seuratObj@meta.data$cci_celltype %>% table() 



##### potential ligands
receiver = c("Post_Endothelial")

expressed_genes_receiver <- get_expressed_genes(receiver, seuratObj, pct = 0.10)

all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)

potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()


sender_celltypes <- c("Sen_Malignant")
list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratObj, 0.10)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()
potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 


length(expressed_genes_sender)
length(potential_ligands)
length(potential_ligands_focused)


#### RUN NICHENETR
ligand_target_matrix = qread('Fig5/data/5b/ligand_target_matrix_endo.qs')
weighted_networks = qread('Fig5/data/5b/weighted_networks_postendo.qs')

geneset = c('CDH2', 'CDH5', 'OCLN', 'CDH13', 'RAMP2', 'NRG3', 'GJB6', 'IL1RAPL1', 'LRRC15', 'CBLN4',
            'PRICKLE1', 'CTTN', 'DKK1', 'FLRT2', 'NEGR1', 'CNTNAP2', 'IL1RAPL2', 'GREM1', 'GJA5', 'GJB2',
            'GPM6A', 'CDH20', 'AMIGO2', 'ITGA6', 'KDR', 'LSAMP', 'CAPRIN1', 'MAP1B', 'CLDN11', 'PLXNA2',
            'WNT4', 'CRTAC1', 'PARD3', 'NLGN4X', 'ROBO2', 'SFRP1', 'DST', 'TEK', 'THBS1', 'KIRREL3', 
            'ACTN1', 'CLDN1', 'ESAM', 'NTN1', 'CADM1', 'ITGB1', 'PTPRF', 'CD80', 'CD34', 'MMP9', 'AQP3')

background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

pot_genes = unique(potential_ligands)
pot_genes = intersect(unique(weighted_networks$lr_sig$from), pot_genes)
length(pot_genes) ; length(potential_ligands)


ligand_activities = predict_ligand_activities(geneset = geneset, 
                                              background_expressed_genes = background_expressed_genes,
                                              ligand_target_matrix = ligand_target_matrix, 
                                              potential_ligands = pot_genes)
ligand_activities = ligand_activities %>% arrange(-auroc)
ligand_activities

qsave(ligand_activities, 'Fig5/data/5b/ligand_activites_endo.qs', preset='fast')
# ligand_activities = qread('Fig5/data/5b/ligand_activites_endo.qs')


## -------------------
## NICHENETR PLOTS
## ------------------
sorted_mat = ligand_target_matrix[intersect(rownames(ligand_target_matrix), rownames(seuratObj@assays$RNA$scale.data)),]

rna_data = seuratObj@assays$RNA$scale.data
annot_data = data.frame(colnames(rna_data),seuratObj$cci_celltype)

Sen_Malignant = data.frame(t(rna_data[,which(seuratObj$cci_celltype == 'Sen_Malignant')]))
Sen_Malignant = apply(Sen_Malignant,2,mean)
Sen_Endothelial = data.frame(t(rna_data[,which(seuratObj$cci_celltype == 'Sen_Endothelial')]))
Sen_Endothelial = apply(Sen_Endothelial,2,mean)
Malignant = data.frame(t(rna_data[,which(seuratObj$cci_celltype == 'malignant')]))
Malignant = apply(Malignant,2,mean)
Non_Malignant = data.frame(t(rna_data[,which(seuratObj$cci_celltype == 'Non_Malignant')]))
Non_Malignant = apply(Non_Malignant,2,mean)
Post_Endothelial = data.frame(t(rna_data[,which(seuratObj$cci_celltype == 'Post_Endothelial')]))
Post_Endothelial = apply(Post_Endothelial,2,mean)
On_Endothelial = data.frame(t(rna_data[,which(seuratObj$cci_celltype == 'On_Endothelial')]))
On_Endothelial = apply(On_Endothelial,2,mean)
Pre_Endothelial = data.frame(t(rna_data[,which(seuratObj$cci_celltype == 'Pre_Endothelial')]))
Pre_Endothelial = apply(Pre_Endothelial,2,mean)

color = colorRampPalette(rev(brewer.pal(n = 16, name ="RdYlBu")))(100)
expression_df = data.frame(Sen_Malignant,Sen_Endothelial,
                           Post_Endothelial, On_Endothelial, Pre_Endothelial,
                           Malignant, Non_Malignant
)

rownames(expression_df) = rownames(seuratObj@assays$RNA$scale.data)

### best ligands
pat.mal = subset(pat, subset = celltype_f %in% c('Endothelial'))
pat.mal = ScaleData(pat.mal, features = geneset)
end.avg.exp <- AverageExpression(pat.mal, group.by = "cci_celltype", features = geneset, slot = "data")
end.avg.exp = end.avg.exp$RNA
end.avg.exp = t(end.avg.exp)
end.avg.exp = t(scale(end.avg.exp))

vis_target_tumor_expression_scaled = t(end.avg.exp)
rownames(vis_target_tumor_expression_scaled) = gsub('-','_', rownames(vis_target_tumor_expression_scaled))

#############

best_upstream_ligands = c('IL1B','VCAM1','CXCL8','CXCL1','FGF2','MMP12','FN1','MMP1','MMP2','TIMP1')

active_ligand_target_links_df = best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,geneset = geneset, ligand_target_matrix = ligand_target_matrix, n = 50000) %>% bind_rows() %>% drop_na()
ctive_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df,
                                                                 ligand_target_matrix = ligand_target_matrix, cutoff = 0.2)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()

rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23


## Fig5b AUPR ligand activity
ligand_aupr_matrix = ligand_activities %>% dplyr::select(aupr) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
vis_ligand_aupr = ligand_aupr_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("aupr")
p_ligand_aupr = vis_ligand_aupr %>% make_heatmap_ggplot("ligands","Ligand activity", color = "orange",
                                                        legend_position = "top", x_axis_position = "top",
                                                        legend_title = "AUPR\n(target gene prediction ability)") +
  theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "white",mid = 'white',high = "red")

p_ligand_aupr

ggsave('Fig5/figures/Fig5b_postendo_ligand_AUPR.pdf', width = 5, height=8, dpi=300)


## Fig5b ligand expression
annot_data2 = data.frame(colnames(rna_data),seuratObj$celltype_f)

expression_df_CAF <- t(rna_data[best_upstream_ligands,]) %>% data.frame() %>%
  tibble::rownames_to_column("cell") %>% as_tibble() %>%
  inner_join(annot_data2 %>% dplyr::rename('cell' = colnames.rna_data., 'cci_celltype' = seuratObj.celltype_f) %>%
               dplyr::select(cell,cci_celltype), by =  "cell")

aggregated_expression_CAF <- expression_df_CAF %>% group_by(cci_celltype) %>%
  dplyr::select(-cell) %>% summarise_all(mean)

aggregated_expression_df_CAF <- aggregated_expression_CAF %>% dplyr::select(-cci_celltype) %>% t() %>%
  magrittr::set_colnames(aggregated_expression_CAF$cci_celltype) %>%
  data.frame() %>% rownames_to_column("ligand") %>% as_tibble()

aggregated_expression_matrix_CAF <- aggregated_expression_df_CAF %>% dplyr::select(-ligand) %>% as.matrix() %>%
  magrittr::set_rownames(aggregated_expression_df_CAF$ligand)


order_tumors <- c('sen_malignant','sen_endothelial','malignant', 'Endothelial')

vis_ligand_tumor_expression <- aggregated_expression_matrix_CAF[rev(best_upstream_ligands), order_tumors]

vis_ligand_tumor_expression_q <- t(
  apply(vis_ligand_tumor_expression, 1, function(x) rescale(x, to = c(0, 1)))
)
p_ligand_tumor_expression = vis_ligand_tumor_expression_q %>%
  make_threecolor_heatmap_ggplot("Ligand", "Cell type", low_color = "#313695",mid_color = "white", mid = 0.5, high_color = "#A50026", legend_position = "top", x_axis_position = "top" , legend_title = "Scaled expression\n(averaged over\nsingle cells)") +
  theme(axis.text.x = element_text(face = "italic"),
        axis.text.y = element_text(size = 15))

p_ligand_tumor_expression

ggsave('Fig5/figures/Fig5b_postendo_ligand+exp.png', width = 5, height=8, dpi=300)
ggsave('Fig5/figures/Fig5b_postendo_ligand+exp.pdf', width = 5, height=8, dpi=300)


## Fig5b ligand-target gene
order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))
order_targets = setdiff(order_targets, c('ACTA2','MMP9','RHCG','S100A8','S100A9'))

order_targets2 = order_targets
length(intersect(order_targets, order_targets2)) ; length(order_targets2)


vis_ligand_target <- t(active_ligand_target_links[order_targets2,order_ligands])

p_ligand_target_network = make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                                              color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")

p_ligand_target_network

ggsave('Fig5/figures/Fig5b_postendo_ligand+target+network.png', width = 10, height=5, dpi=300)
ggsave('Fig5/figures/Fig5b_postendo_ligand+target+network.pdf', width = 10, height=5, dpi=300)


## Fig5b target genes heatmap
ordered = intersect(colnames(vis_ligand_target),colnames(vis_target_tumor_expression_scaled))
ordered = setdiff(order_targets2, c('MMP9','RHCG','S100A8','S100A9'))

order_tumors2 = c('Post_Endothelial','On_Endothelial','Pre_Endothelial')

p_target_tumor_scaled_expression <- (
  apply(vis_target_tumor_expression_scaled[order_tumors2,ordered], 2, function(x) rescale(x, to = c(0, 1)))
)

p_target_tumor_scaled_expression  = p_target_tumor_scaled_expression %>%
  make_threecolor_heatmap_ggplot("Cell type","Target", low_color = "#313695",mid_color = 'white', mid = 0.5, high_color = "#A50026", legend_position = "top", x_axis_position = "top" , legend_title = "Scaled expression\n(averaged over\nsingle cells)") + theme(axis.text.x = element_text(face = "italic"))
p_target_tumor_scaled_expression

ggsave('Fig5/figures/Fig5b_postendo_ligand+target+exp.png', width = 10, height=3, dpi=300)
ggsave('Fig5/figures/Fig5b_postendo_ligand+target+exp.pdf', width = 10, height=3, dpi=300)


## Fig5b ligand-Receptor gene
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig)

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both")

# checked expression
sel.no = c('IL1RAP', 'MMP1','EGFR','FLT4','COL13A1','ITGA8','CD44','EZR')
sel.re = setdiff(unique(dp_data$features.plot), sel.no)

vis_ligand_receptor_network <- vis_ligand_receptor_network[rev(sel.re), rev(best_upstream_ligands)]

p_ligand_receptor_network = (make_heatmap_ggplot(t(vis_ligand_receptor_network),
                                                 y_name = "Ligands", x_name = "Receptors",
                                                 color = "mediumvioletred", legend_title = "Prior interaction potential"))
p_ligand_receptor_network

ggsave('Fig5/figures/Fig5b_postendo_ligand+receptor+network.png', width = 8, height=5, dpi=300)
ggsave('Fig5/figures/Fig5b_postendo_ligand+receptor+network.pdf', width = 8, height=5, dpi=300)




### pathway
library(msigdb)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)

msigdb_go <- msigdbr(species = "Homo sapiens", collection = "C5")
msigdb_go <- msigdb_go[grep('^GO', msigdb_go$gs_name),]
go_term2gene <- msigdb_go[, c("gs_name", "gene_symbol")]

go_emap2 = enricher(order_targets,
                    TERM2GENE = go_term2gene,
                    pvalueCutoff = 0.05)
go_emap2@result = go_emap2@result %>% filter(p.adjust <= 0.05)

go_emap2@result$ID[1:100]

sel.p = c(
  "GOBP_REGULATION_OF_CELL_JUNCTION_ASSEMBLY",
  "GOBP_CELL_MATRIX_ADHESION" ,
  "GOBP_FOCAL_ADHESION_ASSEMBLY" ,
  "GOBP_APICAL_JUNCTION_ASSEMBLY",
  "GOBP_TIGHT_JUNCTION_ORGANIZATION",
  "GOMF_INTEGRIN_BINDING",
  "GOBP_REGULATION_OF_VASCULAR_PERMEABILITY" ,
  "GOMF_CADHERIN_BINDING" ,
  "GOCC_ADHERENS_JUNCTION",
  "GOBP_CELL_JUNCTION_DISASSEMBLY"
)

barplot(go_emap2, showCategory = sel.p)


ggsave('Fig5/figures/Fig5b_postendo_target_pathway.png', width=12, height = 8, dpi=300)
ggsave('Fig5/figures/Fig5b_postendo_target_pathway.pdf', width=12, height = 8, dpi=300)

