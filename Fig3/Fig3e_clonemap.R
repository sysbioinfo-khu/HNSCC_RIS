
#### conda activate cm

library(dplyr)
library(cloneMap)

setwd('/data_hdd4/yuseong/cancer_aging/analysis/h.github/Fig3')

cluster_info = read.table('data/3e-f/conipher/output/Trees/clusterInfo.txt', sep='\t',
                          header=T)

cluster_info <- cluster_info %>%
  mutate(
    ccf_scaled = meanCCF / 100,
    ccf_final = pmin(pmax(ccf_scaled, 0), 1)
  )

ccf_1 = cluster_info %>% filter(SAMPLE == 'HyPR.HN_02_1') %>% select(clusterID, ccf_final)
colnames(ccf_1) = c('clones', 'CCF')
ccf_1$clones = as.character(ccf_1$clones)
ccf_1

ccf_2 = cluster_info %>% filter(SAMPLE == 'HyPR.HN_02_3') %>% select(clusterID, ccf_final)
colnames(ccf_2) = c('clones', 'CCF')
ccf_2$clones = as.character(ccf_2$clones)
ccf_2

tree = readRDS('data/3e-f/conipher/output/Trees/HyPR.HN.02.tree.RDS')

tree_1 = matrix(c('1','3',
                  '1','8',
                  '3','7',
                  '7','2',
                  '7','4',
                  '2','9'), ncol=2, byrow=TRUE)
tree_2 = matrix(c('1','3',
                  '3','7',
                  '7','4',
                  '7','6'), ncol=2, byrow=TRUE)

tree_all = matrix(c('1','3',
                    '1','8',
                    '3','7',
                    '7','2',
                    '7','4',
                    '2','9',
                    '7','6'), ncol=2, byrow=TRUE)

clone.names <- sort(unique( c( tree_all[,1], tree_all[,2] ) ))

clone_colours_example <- c( "#AEEEEE", "#1E90FF", "#AEEEEE", "#A2CD5A", 
                            "#FF6347", "#AEEEEE","#27408B", "#1E90FF" ) 

names(clone_colours_example) <- clone.names

pdf_file <- "figures/Fig3e_CloneMap_1.pdf"
pdf(file = pdf_file, width = 10, height = 7)

cloneMap(tree_1, ccf_1, output.Clone.map.obj = TRUE, clone.cols = clone_colours_example)

dev.off()

pdf_file <- "figures/Fig3e_CloneMap_2.pdf"
pdf(file = pdf_file, width = 10, height = 7)

cloneMap(tree_2, ccf_2, output.Clone.map.obj = TRUE, clone.cols = clone_colours_example)

dev.off()


