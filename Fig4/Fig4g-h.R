
setwd('/data_hdd4/yuseong/cancer_aging/analysis/h.github/')

rm(list = ls())


library(qs)        
library(dplyr)       
library(purrr) 
library(Seurat)   
library(tidygraph)    
library(ggraph)       
library(ggnewscale)  
library(ggplot2)       
library(scales)     
library(msigdbr)   
library(clusterProfiler) 
library(org.Hs.eg.db) 


## load malignant data
mal.so = qread('Fig2/data/2e_g_i/patient_mal_seurat_+GSVA+ME_ssgsea_bicor.qs')
mal.so <- NormalizeData(mal.so)

sen_cl = read.csv('Fig2/data/2e_g_i/patient_mal_clusters.csv',
                  row.names=1)
sen_cl$barcode = rownames(sen_cl)
table(sen_cl$sen_cluster)

meta = mal.so@meta.data
meta$barcode = rownames(meta)

meta_m = left_join(meta, sen_cl, by = 'barcode')

mal.so@meta.data = meta_m
rownames(mal.so@meta.data) = mal.so$barcode

colnames(mal.so@meta.data)

mal.so$sen_cluster = factor(mal.so$sen_cluster, levels = c('C0','C1', 'C2', 'C3', 'C4', 'C5','C6', 'C7','C8', 'C9'))
Idents(mal.so) = mal.so$sen_cluster
DimPlot(mal.so)

## load C6-specific GRN
adj = read.csv('Fig4/data/4g-h/hg38_c6_tf.csv')

adj$TF = adj$source
adj$importance = adj$coef_abs
summary(adj$p)

tfs <- c("EGR1", "POU2F2", "FLI1", "NFATC2")

## extract top 50 target genes
top_k = 50
include_pattern = NULL

edge.df <- map_dfr(tfs, function(tf){
  tmp <- adj %>% dplyr::filter(TF == tf) %>% dplyr::arrange(dplyr::desc(importance))
  keep <- seq_len(min(top_k, nrow(tmp)))
  if (!is.null(include_pattern)) {
    keep <- sort(unique(c(keep, grep(include_pattern, tmp$target))))
  }
  tmp[keep, c("TF","target","importance")] %>%
    dplyr::transmute(from = TF, to = target, weight = importance*5, group = tf)
})


## calculate logFC of target genes
markers = FindMarkers(mal.so, ident.1 = 'C6', logfc.threshold = 0)

seurat_obj = mal.so
acts_mat = read.csv('Fig4/data/4b/sc_mal_VIPER_pivot.csv', row.names = 1)
auc_assay_data <- t(acts_mat)

common_cells <- intersect(colnames(seurat_obj), colnames(auc_assay_data))
auc_assay_data <- auc_assay_data[, common_cells]
seurat_obj <- subset(seurat_obj, cells = common_cells)

seurat_obj[["AUC"]] <- CreateAssayObject(data = auc_assay_data)

seurat_obj <- ScaleData(seurat_obj, assay = "AUC")
DefaultAssay(seurat_obj) <- "AUC"

auc_markers <- FindMarkers(
  seurat_obj,
  ident.1 = 'C6',
  only.pos = TRUE,
  logfc.thershold = 0
)

tf.dat <- data.frame(
  id = tfs,
  type = "TF"
)
head(tf.dat)

auc_markers$id = gsub('Regulon-', '',rownames(auc_markers))
tf.dat = left_join(tf.dat, auc_markers[paste0('Regulon-', unique(tf.dat$id)),])
head(tf.dat)

markers$to = rownames(markers)
edge.df = left_join(edge.df, markers)
head(edge.df)

targets.dat <- edge.df %>%
  dplyr::select(id = to, avg_log2FC = avg_log2FC, p_val_adj = p_val_adj) %>%
  mutate(type = "Target") %>%
  filter(!id %in% tfs)

head(targets.dat)

nodes <- bind_rows(tf.dat, targets.dat) %>%
  mutate(
    significance = ifelse(p_val_adj < 0.05, "Significant", "Not Significant"),
    label = ''
  )

edges <- edge.df %>%
  dplyr::select(from, to, weight)

nodes$label = ifelse(nodes$type == 'TF', nodes$id, '')


#### Fig 4g - TF-target gene NETWORK PLOT

set.seed(3)

## SASP atlas
sasps = c('A2M','AEBP1','CCT2','CLEC3B','COL1A1','COL1A2','COL6A3','CXCL5','CXCL8',
          'DCN','ESYT1','FABP5','FBN1','IGFBP3','IGFBP5','LUM','MMP9','NAMPT','PCOLCE',
          'PLAT','RARRES2','S100A4','SPARC','STC1','THY1','TNC','TNFAIP6','TXLNA','TYMS',
          'UTRN','YEATS4')

gr <- tbl_graph(nodes = nodes[!duplicated(nodes),], edges = edges)
gr <- gr %>%
  activate(nodes) %>%
  mutate(
    tval = as.numeric(avg_log2FC),  
    pval = as.numeric(p_val_adj)  
  )
gr <- gr %>%
  activate(nodes) %>%
  mutate(significance = case_when(
    type == 'TF' ~ 'Significant',
    (id %in% sasps) & (significance == 'Significant') ~ "Significant SASPs",
    significance == "Significant" ~ "Significant",
    TRUE ~ "Not Significant"
  ))

ggraph(gr, layout = 'fr') + 
  geom_edge_link(aes(alpha = weight), color = "grey60") +
  geom_node_point(aes(
    filter = type == 'Target',
    fill = avg_log2FC,           
    size = type,         
    color = significance, 
    stroke = significance, 
  ), shape = 21) +       
  
  scale_fill_gradient2(
    low = "blue", mid = "ghostwhite", high = "red", midpoint = 0,
    name = "Expression - logFC",
    limits = c(-8, 8),       
    oob = scales::squish     
  ) +
  scale_size_manual(values = c("Target" = 5, "TF" = 12), guide = "none") +
  
  new_scale_fill() +   
  geom_node_point(aes(
    filter = type == "TF",
    fill = avg_log2FC,           
    size = type,         
    color = significance, 
    stroke = significance, 
  ), shape = 22) +           
  
  scale_fill_gradient2(
    low = "#91A28C", mid = "ghostwhite", high = "#8F2C37", midpoint = 0,
    name = "TF activity - logFC",
    limits = c(-4, 4),        
    oob = scales::squish
  ) +
  
  geom_node_text(aes(label = label), repulse = TRUE, 
                 fontface = "bold", size = 4) +
  
  scale_color_manual(
    values = c("Significant SASPs" = "red", "Significant" = "black", "Not Significant" = "grey"),
    name = "FDR < 0.05"
  ) +
  scale_discrete_manual(aesthetics = "stroke",
                        values = c("Significant SASPs" = 1.5, "Significant" = 1.5, "Not Significant" = 0),
                        guide = "none") +
  
  scale_edge_alpha(range = c(0.2, 1), guide = "none") +
  
  theme_void()

ggsave('Fig4/figures/Fig4g_patient_sc_TF-Network_top50.pdf', width = 12, height=8, dpi=300)



#### Fig4h - gene set enrichment analysis
tg_pathw = read.table('Fig4/data/4g-h/target_genes_pathway_results.tsv', sep='\t', header=T)

plot_data <- tg_pathw %>%
  mutate(Description = gsub('_',' ',Description)) %>%
  mutate(
    minus_log10_pval = -log10(pvalue),
    Description = factor(Description, levels = Description[order(minus_log10_pval)])
  )

ggplot(plot_data, aes(x = Description, y = Count, fill = minus_log10_pval)) +
  geom_col(width = 0.7) + 
  coord_flip() + 
  scale_fill_gradient(low = "#FFD1D1", high = "#FF0000", name = "-log10(Pvalue)") +
  labs(
    x = NULL,           
    y = "Count"
  ) +
  
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10, face = "bold", color = "black"), 
    panel.grid.major.y = element_blank()
  )

ggsave('Fig4/figures/Fig4h_patient_sc_TF-targets+pathwayenrichment_top50.pdf', width = 8, height=6, dpi=300)


