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
library(scales)


### SET CORES TO USE
library(BiocParallel)
library(RhpcBLASctl)
blas_set_num_threads(30)
omp_set_num_threads(30)


#### load data
pat = qread('Fig2/data/2a-2d/patient_all_processed.qs')
DimPlot(pat, group.by='celltype_f', label=T)

fib.so = qread('Fig6/data/6a-c/patient_seurat_fibroblast_subclusters.qs')
fib.meta = fib.so@meta.data %>% dplyr::select(Fibroblast_sub, barcode)

pat.meta = pat@meta.data
pat.meta$barcode = rownames(pat.meta)

meta = left_join(pat.meta, fib.meta)
rownames(meta) = meta$barcode
meta$Fibroblast_sub = as.character(meta$Fibroblast_sub)
meta$fib.ct = ifelse(is.na(meta$Fibroblast_sub),  meta$celltype_f, meta$Fibroblast_sub)
table(meta$fib.ct)

pat@meta.data = meta

DimPlot(pat, group.by = 'fib.ct')

pat[['cci_celltype']] = pat[['fib.ct']]
table(pat$cci_celltype)

fib.so = subset(fib.so, subset = condition == 'Post')

pat = subset(pat, subset = condition == 'Post')
pat = SetIdent(pat, value = "cci_celltype")


seuratObj <- alias_to_symbol_seurat(pat, "human")
seuratObj@meta.data$cci_celltype %>% table() 


#### potential ligands
receiver = c("F2")
expressed_genes_receiver <- get_expressed_genes(receiver, seuratObj, pct = 0.05)

all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)

potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()


sender_celltypes <- c("sen_malignant")
list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratObj, 0.05)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()
potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 


length(expressed_genes_sender)
length(potential_ligands)
length(potential_ligands_focused)


#### RUN NICHENETR
ligand_target_matrix = qread('Fig6/data/6h/ligand_target_matrix_senmal+F2_post.qs')
weighted_networks = qread('Fig6/data/6h/weighted_networks_senmal+F2_post.qs')


geneset = c('MMP1', 'CDKN1A', 'IL6', 'VEGFA', 'PLAU', 'CCND1', 'MMP3', 'IER3', 'SERPINE1', 'BMP2', 
            'CCL5', 'EDN1', 'DUSP5', 'LIF', 'FGF2', 'F3', 'RIPK2', 'CXCL1', 'MT2A', 'FOSL1', 
            'TNFRSF11B', 'CSF2', 'SDC4', 'HSD11B1', 'CD82', 'GADD45B', 'BIRC3', 'RCAN1', 'CCL7', 'CCR7',
            'CSF3', 'INHBA', 'SLC2A1', 'AKR1B1')

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

qsave(ligand_activities, 'Fig6/data/6h/ligand_activites_F2_filtered_post.qs', preset='fast')
# ligand_activities = qread('Fig6/data/6h/ligand_activites_F2_filtered_post.qs')



## ------------------
## NICHENETR PLOTS
## ------------------
sorted_mat = ligand_target_matrix[intersect(rownames(ligand_target_matrix), rownames(fib.so@assays$RNA$scale.data)),]
sorted_mat = sorted_mat[order(rowSums(sorted_mat[, sen_mal]), decreasing = TRUE), ]
sorted_mat[, sen_mal][1:5,]

rna_data = fib.so@assays$RNA$scale.data
fib.so$cci_celltype = fib.so$Fibroblast_sub
annot_data = data.frame(colnames(rna_data),fib.so$cci_celltype)

C1 = data.frame(t(rna_data[,which(fib.so$cci_celltype == 'F0')]))
C1 = apply(C1,2,mean)
C2 = data.frame(t(rna_data[,which(fib.so$cci_celltype == 'F1')]))
C2 = apply(C2,2,mean)
C3 = data.frame(t(rna_data[,which(fib.so$cci_celltype == 'F2')]))
C3 = apply(C3,2,mean)
C4 = data.frame(t(rna_data[,which(fib.so$cci_celltype == 'F3')]))
C4 = apply(C4,2,mean)
C5 = data.frame(t(rna_data[,which(fib.so$cci_celltype == 'F4')]))
C5 = apply(C5,2,mean)

color = colorRampPalette(rev(brewer.pal(n = 16, name ="RdYlBu")))(100)
expression_df = data.frame(C1,C2,C3,C4,C5)
rownames(expression_df) = rownames(fib.so@assays$RNA$scale.data)


vis_target_tumor_expression_scaled = t(expression_df[rownames(sorted_mat)[1:100],])
colnames(vis_target_tumor_expression_scaled) = rownames(sorted_mat)[1:100]
vis_target_tumor_expression_scaled = scale_quantile(vis_target_tumor_expression_scaled)


#############
best_upstream_ligands = c('IL1B','MMP1','MMP2','MMP14','SPARC')

active_ligand_target_links_df = best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,geneset = geneset, ligand_target_matrix = ligand_target_matrix, n = 50000) %>% bind_rows() %>% drop_na()
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df,
                                                                 ligand_target_matrix = ligand_target_matrix, cutoff = 0.2)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
order_targets = order_targets[!grepl('H1.0', order_targets)]

rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

## Fig6h AUPR ligand activity
ligand_aupr_matrix = ligand_activities %>% dplyr::select(aupr) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
vis_ligand_aupr = ligand_aupr_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("aupr")
p_ligand_aupr = vis_ligand_aupr %>% make_heatmap_ggplot("ligands","Ligand activity", color = "orange",
                                                        legend_position = "top", x_axis_position = "top",
                                                        legend_title = "AUPR\n(target gene prediction ability)") +
  theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "white",mid = 'white',high = "red")

p_ligand_aupr

ggsave('Fig6/figures/Fig6h_Fib_ligand_AUPR_post.pdf', width = 5, height=8, dpi=300)


## Fig6h ligand expression
pat.sen = subset(pat, subset = cci_celltype %in% c('sen_malignant','sen_endothelial', 'malignant','Endothelial'))
sen.avg_exp <- AverageExpression(pat.sen, group.by = "cci_celltype", features = sen_mal, slot = "data")
sen.avg_exp = sen.avg_exp$RNA
sen.avg_exp = t(sen.avg_exp)
sen.avg_exp = t(scale(sen.avg_exp))

ordered.sen = c('sen-malignant','sen-endothelial', 'malignant','Endothelial')
sen.avg_exp_q <- t(apply(sen.avg_exp, 1, function(x) rescale(x, to = c(0, 1))))
p_ligand_tumor_expression  = sen.avg_exp_q[order_ligands,ordered.sen]  %>% as.matrix() %>%
  make_threecolor_heatmap_ggplot("Cell type","Target", low_color = "#313695",mid_color = 'white', mid = 0.5, high_color = "#A50026", legend_position = "top", x_axis_position = "top" , legend_title = "Scaled expression\n(averaged over\nsingle cells)") + theme(axis.text.x = element_text(face = "italic"))
p_ligand_tumor_expression

ggsave('Fig6/figures/Fig6h_Fib_ligand+exp_post.pdf', width = 5, height=8, dpi=300)


## Fig6h ligand-target gene
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
p_ligand_target_network = vis_ligand_target %>%
  make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple", legend_position = "top", x_axis_position = "top",
                      legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")
p_ligand_target_network

ggsave('Fig6/figures/Fig6h_Fib_ligand+target+network_post.png', width = 10, height=5, dpi=300)
ggsave('Fig6/figures/Fig6h_Fib_ligand+target+network_post.pdf', width = 10, height=5, dpi=300)



## Fig6h target genes heatmap
ordered = colnames(vis_ligand_target)
order_cafs = c('C1','C2','C3','C4','C5')

vis_target_tumor_expression_scaled = t(expression_df[colnames(vis_ligand_target),])
colnames(vis_target_tumor_expression_scaled) = colnames(vis_ligand_target)
vis_target_tumor_expression_scaled = scale_quantile(vis_target_tumor_expression_scaled)

p_target_tumor_scaled_expression  = vis_target_tumor_expression_scaled[order_cafs,colnames(vis_target_tumor_expression_scaled) %in% ordered]  %>%
  make_threecolor_heatmap_ggplot("Cell type","Target", low_color ="#313695",mid_color = "white", mid = 0.5, high_color = "#A50026", legend_position = "top", x_axis_position = "top" , legend_title = "Scaled expression\n(averaged over\nsingle cells)") + theme(axis.text.x = element_text(face = "italic"))
p_target_tumor_scaled_expression

ggsave('Fig6/figures/Fig6h_Fib_ligand+target+exp_post.png', width = 10, height=3, dpi=300)
ggsave('Fig6/figures/Fig6h_Fib_ligand+target+exp_post.pdf', width = 10, height=3, dpi=300)


## Fig6h ligand-Receptor gene
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)


lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
vis_ligand_receptor_network = t(vis_ligand_receptor_network)

# checked expression
sel.no = c('DDR1', 'F2R', 'ITGAV', 'ITGB3','PECAM1','SDC2','IL1R1','IL1RAP')
sel.re = setdiff(unique(dp_data$features.plot), sel.no)


vis_ligand_receptor_network = vis_ligand_receptor_network[rev(best_upstream_ligands), sel.re]
p_ligand_receptor_network = vis_ligand_receptor_network %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network

ggsave('Fig6/figures/Fig6h_Fib_ligand+receptor+network_post.png', width = 8, height=5, dpi=300)
ggsave('Fig6/figures/Fig6h_Fib_ligand+receptor+network_post.pdf', width = 8, height=5, dpi=300)


### Fig6h pathway of target genes
library(msigdb)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)

msigdb_all = msigdbr(species = "Homo sapiens")
msigdb_sel = msigdb_all %>% filter(gs_collection == 'H' | grepl('CP:REACTOME|GO',gs_subcollection))
unique(msigdb_sel$gs_subcollection)
unique(msigdb_sel$gs_collection)

db_term2gene <- msigdb_sel[, c("gs_name", "gene_symbol")]

inter_emap = enricher(order_targets,
                      TERM2GENE = db_term2gene,
                      pvalueCutoff = 0.05)
inter_emap@result = inter_emap@result %>% filter(pvalue <= 0.05)

sel.id = c('GOMF_CYTOKINE_ACTIVITY',
           'REACTOME_SIGNALING_BY_INTERLEUKINS',
           'GOMF_CYTOKINE_RECEPTOR_BINDING',
           "HALLMARK_INFLAMMATORY_RESPONSE",
           'GOMF_CXCR_CHEMOKINE_RECEPTOR_BINDING',
           "HALLMARK_IL2_STAT5_SIGNALING" ,
           "REACTOME_CHEMOKINE_RECEPTORS_BIND_CHEMOKINES" ,
           "GOMF_CHEMOKINE_ACTIVITY"  )

barplot(inter_emap, showCategory = sel.id)

ggsave('Fig6/figures/Fig6h_Fib_target_pathway_post.pdf', width=8, height = 5, dpi=300)


