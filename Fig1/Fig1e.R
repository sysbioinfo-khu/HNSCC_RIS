
setwd('/data_hdd4/yuseong/cancer_aging/analysis/h.github/Fig1/')

rm(list=ls())


library(qs)              
library(DESeq2)        
library(ggplot2)       
library(dplyr)          
library(EnhancedVolcano)
library(msigdbr)        
library(clusterProfiler)

#### Extended Data Fig 1a - PCA plot of our celllines 
pca_dat = qread('data/1e/KHUcellline_PCA_data.qs')

a = plotPCA(pca_dat, intgroup="group")+ 
  scale_color_manual(values = c('WT' = 'dodgerblue', 'RT' = 'brown', 'PL' = 'blue', 'SEN' = 'tomato')) +
  theme_bw()
a

ggsave('figures/ExtendedFig1a_KHUcellline_PCAplot.pdf', a, dpi=300, width=8,height=6)




#### Extended Data Fig 1b - DEGs between UR and SEN in our celllines

res_sen_df = read.csv('data/1e/KHUcellline_DEG_SEN_vs_WT.csv', row.names=1)

keyvals <- rep('grey', nrow(res_sen_df))
names(keyvals) <- rep('NS', nrow(res_sen_df))

keyvals[which(res_sen_df$log2FoldChange > 1 & res_sen_df$pvalue < 0.05)] <- 'tomato'
names(keyvals)[which(res_sen_df$log2FoldChange > 1 & res_sen_df$pvalue < 0.05)] <- 'SEN_Upregulated'

keyvals[which(res_sen_df$log2FoldChange < -1 & res_sen_df$pvalue < 0.05)] <- 'royalblue'
names(keyvals)[which(res_sen_df$log2FoldChange < -1 & res_sen_df$pvalue < 0.05)] <- 'WT_Upregulated'


v = EnhancedVolcano(res_sen_df,
                    lab = rownames(res_sen_df),
                    x = 'log2FoldChange',
                    y = 'pvalue',
                    pCutoff = 0.05,
                    FCcutoff = 1,
                    colCustom = keyvals) 
v =v + labs(title = '', subtitle = '')
v

ggsave('figures/ExtendedFig1b_KHUcellline_DEG_volcanoplot.pdf', v, dpi=300, width=8,height=8)




#### Fig 1e - GSEA result between UR and SEN in our celllines
sel.sen_gse = read.table('data/1e/KHUcellline_GSEA_SEN_vs_WT.tsv', sep='\t')

sel.sen = c(
  ## CYTOKINE
  "REACTOME_INTERLEUKIN_10_SIGNALING",    
  "GOBP_ACUTE_INFLAMMATORY_RESPONSE" ,
  "HALLMARK_INFLAMMATORY_RESPONSE",   
  "GOBP_INTERLEUKIN_6_PRODUCTION" ,
  "GOBP_CYTOKINE_MEDIATED_SIGNALING_PATHWAY" ,
  
  ## EMT    
  "GOBP_EXTRACELLULAR_MATRIX_DISASSEMBLY",
  "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX",  
  "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", 
  "GOBP_REGULATION_OF_CELL_CELL_ADHESION",  
  
  ## RESPONSE TO UV
  "GOBP_RESPONSE_TO_UV_A" ,
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_HYPOXIA",
  "GOBP_APOPTOTIC_CELL_CLEARANCE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  
  ## CELL CYCLE    
  "REACTOME_SYNTHESIS_OF_DNA",
  "GOBP_DNA_REPLICATION_INITIATION",  
  "REACTOME_DNA_STRAND_ELONGATION",
  "GOBP_REGULATION_OF_DNA_TEMPLATED_DNA_REPLICATION",    
  "GOBP_CELL_CYCLE_DNA_REPLICATION"
)

sel.sen_gse$Description = factor(sel.sen_gse$Description, levels = rev(sel.sen))

A = ggplot(sel.sen_gse, aes(x = GeneRatio, y = Description)) +
  geom_point(aes(size = count, fill = pvalue), shape = 21, color = "black", stroke = 0.5) +
  facet_grid(. ~ (NES < 0), labeller = as_labeller(c("FALSE" = "Activated", "TRUE" = "Suppressed"))) +
  scale_fill_gradient(high = "#CC3333", low = "#466f9d", trans = "reverse") +
  theme_bw() 

A

ggsave('figures/Fig1e_KHUcellline_Pathway_SEN.pdf', A, dpi=300, width=8,height=8)




