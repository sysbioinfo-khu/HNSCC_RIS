setwd("/data_hdd4/yuseong/cancer_aging/analysis/h.github/")

library(ape)
library(phylogram)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SpatialInferCNV)
library(tidyverse)
library(ggtree)
library(dendextend)
library(qs)
library(reshape2)
library(RColorBrewer)
# library(yarrr)
library(scales)

rm(list=ls())

#===========================================
## LOAD SINGLECELLOBJECT
#===========================================

mal = qread('Fig3/data/3c/patient_epithelial_subset.qs')
mal2 = subset(mal, subset = grep('HyPR-HN-01',Sample))
unique(mal2$Sample)


#===========================================
## LOAD INFERCNV OBJECT AND DENDROGRAM
#===========================================

input_path <- "Fig3/data/3c/infercnv/HyPR-HN-01/" 
ifobj = readRDS(file.path(input_path, 'run.final.infercnv_obj'))
dendrogram = ifobj@tumor_subclusters$hc$all_observations

mytree <- as.phylo(dendrogram)
my.subtrees <- subtrees(mytree)

#Generate a phylogenetic tree with node IDs for downstream manual clone assignment
png(paste0("Fig3/figures/hy01_phylo.png"), width = 4000, height = 10000, res = 300)
plot(mytree, show.tip.label = FALSE)
nodelabels(text = 1:mytree$Nnode, node= 1:mytree$Nnode + Ntip(mytree))
dev.off()


#===================================================================
## CLONE PLOTS
#===================================================================
Node1 <- rbind(SelectingSubTreeData(my.subtrees, 2), SelectingSubTreeData(my.subtrees, 4))
Node2 <- SelectingSubTreeData(my.subtrees, 6)
Node3 <- SelectingSubTreeData(my.subtrees, 8)
Node4 <- SelectingSubTreeData(my.subtrees, 9)

#Merging All Nodes Together
MergedNodes <- rbind(Node1, Node2, Node3, Node4)

#Then drop the Histology column
colnames(MergedNodes) = c("Barcode", "Node")

length(unique(MergedNodes$Barcode))
length(MergedNodes$Barcode)

table(MergedNodes$Node)

#check groupings
grouping = read.table(file.path(input_path,'infercnv_subclusters.observation_groupings.txt'))

table(grouping$Annotation.Group)

#This then renames "Node" to "Clone"
MergedNodes$Node <- ifelse(MergedNodes$Node == "Node_6" , "Clone C",
                           ifelse(MergedNodes$Node == "Node_8" , "Clone A",
                                  ifelse(MergedNodes$Node == "Node_9" , "Clone B",
                                         ifelse(MergedNodes$Node %in% c("Node_2", "Node_4"), "Clone D",
                                                MergedNodes$Node))))
table(MergedNodes$Node)

#Include clone information in the metadata.
mal2@meta.data$Barcode = rownames(mal2@meta.data)

meta = mal2@meta.data
meta_m = left_join(meta, MergedNodes, by = 'Barcode') %>% dplyr::rename('Clone' = Node)
rownames(meta_m) = meta_m$Barcode
mal2@meta.data = meta_m
mal2@meta.data$Clone = ifelse(is.na(mal2$Clone), 'NT', mal2$Clone)

DimPlot(mal2, group.by = 'Clone', label=T)

write.csv(mal2@meta.data[c('Clone')], 'Fig3/data/3c/clone_analysis/hy01_clone_metadata.csv')

#--------------------------------------------------

#Information update
Idents(mal2) <- factor(mal2@meta.data$Clone, levels = sort(names(table(mal2@meta.data$Clone))))

#Define a color panel for tumor clone
clone.color <- c('#F48B9A', '#5CA484', '#CFB177', '#C9C3BD', 'white')
names(clone.color) <- c(names(table(mal2$Clone)))

#Define a color panel for condition
site.color = c('grey','springgreen4', 'orange')
names(site.color) = c('Pre','On','Post')

#Clustering & flip base on the clone & Manually verify the definition of the tumor clone.
allnode = mytree$tip.label
cluster = MergedNodes$Node
names(cluster) = MergedNodes$Barcode

groupInfo <- lapply(levels(as.factor(cluster)), function(cc, cluster, allnode){
  return(names(cluster)[cluster == cc])
},	cluster, allnode)					 
names(groupInfo) <- names(table(cluster))

MPT = groupOTU(mytree, groupInfo)
gp1 <- ggtree(MPT, aes(color = group)) + theme_tree() + geom_tippoint(alpha = 0.7)
gp1 <- gp1 + scale_color_manual(values=clone.color[grep("Clone", names(clone.color))])
gp1 + geom_text(aes(label = node))

#Calculate the tissue composition for each subclone
row_anno_clean <- mal2@meta.data[,c("condition", "Clone")]
row_anno_clean$condition <- factor(row_anno_clean$condition, levels = c("Pre", "On", "Post"))
row_anno_clean <- row_anno_clean[grep("Clone", row_anno_clean$Clone),]
clone_composition <- table(row_anno_clean$condition, row_anno_clean$Clone)
clone_composition <- t(clone_composition[,grep("Clone", colnames(clone_composition))])

clone_composition_df <- data.frame(matrix(as.vector(clone_composition), ncol = ncol(clone_composition)))
rownames(clone_composition_df) <- rownames(clone_composition)
colnames(clone_composition_df) <- colnames(clone_composition)
gp1 + geom_text(aes(label=node), hjust=0)

#Manually specify the parent node of each tumor clone to add the pie chart.
clone_composition_df$node <-  c(148, 151, 153,  154)
pie_df <- nodepie(clone_composition_df, cols=1:3, alpha=0.9, color = site.color ) 

gp_combo1 <- inset(gp1, pie_df, height=0.1, width = 0.1)
gp_combo2 <- inset(gp1, pie_df, height=0.1, width = 0.1, hjust = -0.035, vjust = 0.03)

dd = gp_combo1@data[which(gp_combo1@data$isTip),]

#Get the order of spot for downstream heatmap
spot_order1 <- dd$label[order(dd$y, decreasing=TRUE)]

#Load the infercnv output
pre.obj = readRDS(file.path(input_path, 'preliminary.infercnv_obj'))
data1 = pre.obj@expr.data
data1 = data1[,WhichCells(mal2, expression = grep('Clone', Clone))]
cnv <- t(data1)
data <- cnv[MergedNodes[grep("Clone", (MergedNodes$Node)),]$Barcode,]

#Remove the genes on chromosomes X and Y.
col_anno <- read.delim("Fig3/data/3c/infercnv/gencode.v48.primary_assembly.annotation.gtf.gen_pos.txt", head = F)
col_anno <- col_anno[col_anno$V1 %in% colnames(cnv),]
col_anno <- col_anno[col_anno$V2 %in% paste0("chr", 1:22),]
rownames(col_anno) <- col_anno[,1]

col_anno_clean <- data.frame(col_anno[,1:2], row.names = rownames(col_anno))
intersect_gene <- intersect(colnames(cnv), rownames(col_anno_clean))
cnv_hdata_clean <- cnv[, intersect_gene]
col_anno_clean <- col_anno_clean[rownames(col_anno_clean) %in% intersect_gene,]
col_anno_clean <- data.frame(Chromosome = col_anno_clean[,2], row.names = rownames(col_anno_clean))
col_anno_clean$Chromosome <- factor(col_anno_clean$Chromosome, levels = paste0("chr", c(1:22)))

tb = table(col_anno_clean);
tb = tb[ paste0("chr", 1:22)]
seg.chr = cumsum(tb)

#Create color panel for heatmap
library(pheatmap)
library(pals) 
pals::cols25()[1:22] -> chr.color
names(chr.color) = paste0("chr", c(1:22))

aCol = list(Chromosome = chr.color, site = site.color, Clone = clone.color)
data_clean <- data[,rownames(col_anno_clean)]

#Generate CNV heatmap
mid_chr <- cumsum(tb)			
mid_tb <- round(tb/2)			
mid_chr_mid <- round(mid_chr/2)
mid_chr_mid[1] <- round(mid_chr[1]/2)
mid_chr_mid[2:length(mid_chr_mid)] <- round(mid_chr[-1] - mid_tb[-1])

data_clean_map <- data_clean

for (k in 1:length(colnames(data_clean_map))){
  colnames(data_clean_map)[k] <- " "
}
colnames(data_clean_map)[mid_chr_mid] <- gsub("chr", "", names(seg.chr))

pheatmap::pheatmap(data_clean_map[spot_order1,], scale = 'none', color = colorRampPalette(c("navy", "white", "firebrick"))(50), angle_col = "90", 
                   cluster_cols = F, gaps_row = NULL, gaps_col = seg.chr, annotation_colors = aCol, fontsize = 8,
                   annotation_row = row_anno_clean,
                   annotation_col = NA, 
                   cluster_rows = F, show_colnames = T, show_rownames = F, border_color = NA,
                   width = 10,
                   height = 5,
                   cutree_col = 2,
                   file = paste0("Fig3/figures/hy01_Figure 2_pt1_heatmap_raw.png") 	)

#Merge the CNV from spot to clone level
dim(data)
dim(MergedNodes)

MergedNodes_clean <- MergedNodes[ MergedNodes$Barcode %in% rownames(data), ]
rownames(MergedNodes_clean) <- MergedNodes_clean$Barcode
length(rownames(data) == rownames(MergedNodes_clean))

data_merge <- data.frame(matrix(1, nrow = length(table(MergedNodes_clean$Node)), ncol = ncol(data)))
rownames(data_merge) <- names(table(MergedNodes_clean$Node))
colnames(data_merge) <- colnames(data)

for (i in 1:ncol(data_merge)){
  data_merge[,i] <- tapply(data[,i], MergedNodes_clean$Node, mean)
}
data_merge <- data_merge[,rownames(col_anno_clean)]

fig_thres <- c(max((data_merge - 1)), min((data_merge - 1)))
fig_thres <- min(abs(fig_thres))
data_merge[data_merge > 1 + fig_thres] <- 1 + fig_thres
data_merge[data_merge < 1 - fig_thres] <- 1 - fig_thres

#Load the infercnv output of reference
data2 = pre.obj@expr.data
refs = names(mal2$Clone)[is.na(mal2$Clone)]
data2 = data2[,refs]
dim(data2)
cnv_ref <- t(data2)

#Alternatively, use a clean reference without any CNV. This node should be clean, as the root node needs to be specified for subsequent evolutionary analysis.
cnv_ref_clean <- matrix(1, nrow = nrow(data_merge), ncol = ncol(data_merge))
rownames(cnv_ref_clean) <- rownames(data_merge)
colnames(cnv_ref_clean) <- colnames(data_merge)
data_merge_ref <- rbind(data_merge, "NT" = apply(cnv_ref_clean, 2, mean))

#Remove those week CNV signal, here we used 1.02 at cutoff
data_merge_ref_clean <- data_merge_ref
data_merge_ref_clean[abs(log(data_merge_ref_clean, 2)) < log(1.02, 2)] <- 1
data_dist <- dist(data_merge_ref_clean)

#NJ tree construction at subclone level
library(ggtree)
library(phangorn)
treeNJ2  <- NJ(data_dist)
my_tree <- as.phylo(treeNJ2)
allnode = my_tree$tip.label
cluster = allnode
names(cluster) = allnode

groupInfo <- lapply(levels(as.factor(cluster)), function(cc, cluster, allnode){
  return(names(cluster)[cluster == cc])
},	cluster, allnode)					 
names(groupInfo) <- names(table(cluster))

#Specify the root clone as NT.
treeNJ2 <- (root(treeNJ2, which(treeNJ2$tip.label == "NT")))
MPT5 = groupOTU(treeNJ2, groupInfo)

#Calculate the spot number
dot_size <- 3 * (log(as.numeric(table(mal2@meta.data$Clone)[grep("Clone|NT", names(table(mal2@meta.data$Clone)))]))) - 3
dot_size[dot_size > 18] <- 18
dot_size <- dot_size + 8

#Replace the NT color
clone.color["NT"] <- "darkgrey"
gp5 <- ggtree(MPT5, aes(color = group), size = 2) + theme_tree() #+ geom_tippoint(alpha = 0.7)
gp5 <- gp5 + scale_color_manual(values=clone.color[grep("Clone|NT", names(clone.color))])
gp5 <- gp5 + geom_tippoint(size = dot_size)
gp5$data$label <- gsub("Clone", "", gp5$data$label)
gp5 <- gp5 + geom_tiplab(hjust = 0.55, vjust = 0.45, color = "black", size = 14) +  theme(legend.position="none")

#Add NT spot number
clone_composition_df_ref <- rbind(clone_composition_df, "NT" = 0)
temp <- table(mal2@meta.data$Clone, mal2@meta.data$condition)[c("NT"),]
clone_composition_df_ref["NT", match(names(temp), colnames(clone_composition_df_ref))] <- temp

gp1 + geom_text(aes(label=node), hjust=0)					 
clone_composition_df_ref$node <- c(1:nrow(clone_composition_df_ref))
pie_df <- nodepie(clone_composition_df_ref, cols=1:3, alpha=0.9, color = site.color ) 

gp_combo5 <- inset(gp5, pie_df, height=0.12, width = 0.2,  hjust = -0.2, vjust = 0.05) + xlim(NA, max(gp5$data$x) + 0.2)
dis_node <- gp5$data[gp5$data$isTip,]
gp_combo5a <- inset(gp5, pie_df, height=0.12, width = 0.2, hjust = c(-0.6 - max(dis_node$x) + dis_node$x),  vjust = 0.05) +  xlim(NA, max(dis_node$x) + 0.6)

pdf("Fig3/figures/FigS3b_hy01_Figure 2_pt1_part1.pdf", 6, 8)
print(gp_combo5a)
dev.off()


#CNV heatmap at subclone level
d5 = fortify(MPT5)
dd = subset(d5, isTip)
spot_order <- dd$label[order(dd$y, decreasing=TRUE)]
spot_order <- as.character(gp5$data$group[gp5$data$isTip == TRUE])[order(gp5$data$y[gp5$data$isTip == TRUE], decreasing = T)]

row_anno_clean_merge <- clone_composition_df[,1:2]
colnames(row_anno_clean_merge)[1:2] <- "Clone"
row_anno_clean_merge[,1:2] <- rownames(row_anno_clean_merge)
row_anno_clean_merge <- as.data.frame(row_anno_clean_merge[,1])
colnames(row_anno_clean_merge) <- "Clone"
rownames(row_anno_clean_merge) <- row_anno_clean_merge[,1]
row_anno_clean_merge[nrow(row_anno_clean_merge)+1,] <- "NT"
rownames(row_anno_clean_merge)[nrow(row_anno_clean_merge)] <- "NT"
data_merge <- data_merge[,rownames(col_anno_clean)]

data_merge_map <- data_merge_ref
colnames(data_merge_map)[1:length(data_merge_map)] <- " "
colnames(data_merge_map)[mid_chr_mid] <- gsub("chr", "", names(seg.chr))
data_merge_map[abs(log(data_merge_map, 2)) < log(1.02, 2)] <- 1

aCol2 <- aCol
aCol2$Clone[names(aCol2$Clone) == "NT"] <- "darkgrey"
aCol2$Clone <- aCol2$Clone[grep("Clone|NT", names(aCol2$Clone))]

pheatmap::pheatmap(data_merge_map[spot_order,], scale = 'none', color = colorRampPalette(c("navy", "white", "firebrick"))(50), angle_col = "90", 
                   cluster_cols = F, gaps_row = NULL, gaps_col = seg.chr, annotation_colors = aCol2, fontsize = 8,
                   annotation_row = row_anno_clean_merge,
                   annotation_col = NA, 
                   cluster_rows = F, show_colnames = T, show_rownames = F, border_color = NA,
                   width = 6,
                   height = 3.5,
                   cutree_col = 2, 
                   file = "Fig3/figures/FigS3b_hy01_Figure 2_pt1_part2.png")

# expand range of 0
data_processed <- data_merge_map
data_processed[data_processed >= 0.96 & data_processed <= 1.04] <- 1.0

min_range <- 0.7
max_range <- 1.3

custom_legend_breaks <- c(0.8, 1.0, 1.2)
custom_legend_labels <- c("0.8", "1.0", "1.2")

custom_breaks_range <- seq(min_range, max_range, length.out = 50 + 1)
custom_colors <- colorRampPalette(c("navy", "white", "firebrick"))(50)


pheatmap::pheatmap(
  data_processed[spot_order,],
  scale = 'none',
  color = custom_colors,
  breaks = custom_breaks_range,            
  legend_breaks = custom_legend_breaks,     
  legend_labels = custom_legend_labels,
  angle_col = "90",
  cluster_cols = F,
  gaps_row = NULL,
  gaps_col = seg.chr,
  annotation_colors = aCol2,
  fontsize = 8,
  annotation_row = row_anno_clean_merge,
  annotation_col = NA,
  cluster_rows = F,
  show_colnames = T,
  show_rownames = F,
  border_color = NA,
  width = 6,
  height = 3.5,
  cutree_col = 2,
  file = "Fig3/figures/FigS3b_hy01_Figure 2_pt1_part2_expandedWhite.png"
)

########################
data_merge_map <- data_merge_ref
# colnames(data_merge_map)[1:length(data_merge_map)] <- " "
colnames(data_merge_map)[mid_chr_mid] <- gsub("chr", "", names(seg.chr))
data_merge_map[abs(log(data_merge_map, 2)) < log(1.04, 2)] <- 1
data_merge_map[1:4, 1:4]

cl.gns = list()

for (gn in colnames(data_merge_map)) {
  for (cl in rownames(data_merge_map)) {
    if (data_merge_map[cl,gn] > 1) {
      tmp = cl.gns[[cl]]$gain
      cl.gns[[cl]]$gain = c(tmp, gn)
    } else if (data_merge_map[cl,gn] < 1) {
      tmp = cl.gns[[cl]]$loss
      cl.gns[[cl]]$loss = c(tmp, gn)
    } else { 
      next
    } 
  }
}

###

library(sigminer)
data(cytobands.hg38)
cytobands.hg38[1:3,]

cytob = cytobands.hg38 %>% dplyr::rename('chr' = chrom)
head(cytob)

gene_info <- read.delim("Fig3/data/3c/infercnv/gencode.v48.primary_assembly.annotation.gtf.gen_pos.txt", head = F)
gene_info <- gene_info[gene_info$V1 %in% colnames(cnv),]
gene_info <- gene_info[gene_info$V2 %in% paste0("chr", 1:22),]
rownames(gene_info) <- gene_info[,1]
colnames(gene_info) = c('gene', 'chr','start', 'end')

gene_info$band = ''
for (i in seq(1:nrow(gene_info))) {
  chr.hmm = gene_info[i, 'chr']
  start.hmm = gene_info[i,'start']
  end.hmm = gene_info[i,'end']
  
  tp1 = cytob[(cytob$chr == chr.hmm) & (cytob$start <= start.hmm) & (cytob$end >= end.hmm),]$band
  tp2 = cytob[(cytob$chr == chr.hmm) & (cytob$start <= end.hmm) & (cytob$end >= end.hmm),]$band
  tp3 = cytob[(cytob$chr == chr.hmm) & (cytob$start <= start.hmm) & (cytob$end >= start.hmm),]$band
  tp4 = cytob[(cytob$chr == chr.hmm) & (cytob$start >= start.hmm) & (cytob$end >= start.hmm) 
              & (cytob$start <= end.hmm)  & (cytob$end <= end.hmm),]$band
  
  bands = unique(c(tp1, tp2, tp3, tp4))
  bands = sapply(bands, function(x) {paste0(gsub('chr', '', chr.hmm), x)})
  
  gene_info[i, 'band'] = paste0(bands, collapse='|') 
  
}

cl.gns.df = data.frame('Clone' = '',
                       'gene' = '',
                       'state' = '',
                       'band' = '')

n=1
for (cl.n in sort(names(cl.gns))) {
  tmp.gns.g = cl.gns[[cl.n]]$gain
  tmp.gns.l = cl.gns[[cl.n]]$loss
  for (gn in tmp.gns.g) {
    if (gn %in% gene_info$gene) {
      cl.gns.df[n,'Clone'] = cl.n
      cl.gns.df[n,'gene'] = gn
      cl.gns.df[n,'state'] = 'gain'
      cl.gns.df[n, 'band'] = gene_info[gene_info$gene == gn, 'band']
      n = n+1
    }
  }
  for (gn in tmp.gns.l) {
    if (gn %in% gene_info$gene) {
      cl.gns.df[n,'Clone'] = cl.n
      cl.gns.df[n,'gene'] = gn
      cl.gns.df[n,'state'] = 'loss'
      cl.gns.df[n, 'band'] = gene_info[gene_info$gene == gn, 'band']
      n = n+1
    }
  }
}

write.csv(cl.gns.df, 'Fig3/data/3c/clone_analysis/hy01_clone_genes_band_cnv.csv', row.names=F)

cl.gns.df = read.csv('Fig3/data/3c/clone_analysis/hy01_clone_genes_band_cnv.csv')
cl.b.df = cl.gns.df %>%
  separate_rows(band, sep = "\\|") %>%
  group_by(Clone,state,band) %>%
  summarise(cnt=n())
write.csv(cl.b.df, 'Fig3/data/3c/clone_analysis/hy01_clone_band_cnv.csv', row.names=F)

