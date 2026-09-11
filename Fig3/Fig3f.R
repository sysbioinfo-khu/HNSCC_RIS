
library(MutationalPatterns)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(dplyr)
library(igraph)

rm(list=ls())

setwd('/data_hdd4/yuseong/cancer_aging/analysis/h.github/Fig3')


##############################
### single cell clones
##############################

rm(list=ls())

ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)

#### mk clone vcf
get_advanced_lineages <- function(graph) {
  all_nodes <- V(graph)$name
  out_degrees <- igraph::degree(graph, mode = "out")
  
  branch_nodes <- names(out_degrees[out_degrees > 1])
  
  leaves <- names(out_degrees[out_degrees == 0])
  
  lineages <- list()
  
  for (leaf in leaves) {
    path_nodes <- c(leaf)
    current <- leaf
    
    while(TRUE) {
      parents <- neighbors(graph, current, mode = "in")$name
      if (length(parents) == 0) break
      
      parent <- parents[1]
      if (igraph::degree(graph, parent, mode = "out") == 1) {
        path_nodes <- c(path_nodes, parent)
        current <- parent
      } else {
        break
      }
    }
    lineages[[paste0("Lineage_Node_", leaf)]] <- as.numeric(path_nodes)
  }
  return(lineages)
}


## single cell
tree_obj = readRDS('Fig3/data/3e-f/conipher/output/Trees/HyPR.HN.02.tree.RDS')
tree_table <- read.table("Fig3/data/3e-f/conipher/output/Trees/treeTable.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

edges = tree_obj$graph_pyclone$default_tree

g <- graph_from_data_frame(edges, directed = TRUE)

merged_clones <- get_advanced_lineages(g)
merged_clones


merge_lineages <- function(lineages, groups) {
  used <- unlist(groups, use.names = FALSE)
  stopifnot(all(used %in% names(lineages)))
  out <- lineages[setdiff(names(lineages), used)]
  for (nm in names(groups))
    out[[nm]] <- sort(unique(unlist(lineages[groups[[nm]]])))
  out
}

merged_clones <- merge_lineages(
  merged_clones,
  list(PreOnly = c("Lineage_Node_8", "Lineage_Node_9"))
)
print(merged_clones)


vcf_list <- list()
for (lineage_name in names(merged_clones)) {
  cluster_ids <- merged_clones[[lineage_name]]
  
  merged_snvs <- tree_table %>%
    filter(treeCLUSTER %in% cluster_ids,
           nchar(REF) == 1, nchar(ALT) == 1)
  
  cat(lineage_name, ": clusters", paste(cluster_ids, collapse = ","),
      "->", nrow(merged_snvs), "SNVs\n")
  
  if (nrow(merged_snvs) >= 10) {
    gr <- GRanges(seqnames = paste0("chr", merged_snvs$CHR),
                  ranges   = IRanges(merged_snvs$POS, merged_snvs$POS),
                  ref      = merged_snvs$REF,
                  alt      = merged_snvs$ALT)
    seqlevelsStyle(gr) <- "UCSC"
    genome(gr) <- "hg38"                       # ★ hg19 아님
    gr <- keepStandardChromosomes(gr, pruning.mode = "tidy")
    gr <- gr[seqnames(gr) %in% paste0("chr", c(1:22, "X", "Y"))]
    
    vcf_list[[paste0("HyPR.HN.02_", lineage_name)]] <- gr
  }
}

stopifnot(length(vcf_list) > 0)
mut_mat_merged <- mut_matrix(vcf_list  = vcf_list,
                             ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
colSums(mut_mat_merged)     

head(mut_mat_merged)


signatures <- get_known_signatures(genome = 'GRCh38')

fit_res <- fit_to_signatures(mut_mat_merged, signatures)

contributions_bootstrapped <- fit_to_signatures_bootstrapped(
  mut_mat_merged, 
  signatures, 
  n_boots = 100, 
  method = "strict"
)

plot_bootstrapped_contribution(contributions_bootstrapped, 
                               mode = "relative", 
                               plot_type = "dotplot") + coord_flip()

ggsave('figures/Fig3f_WES_SBSplt_clones_strict_pre_merged.pdf',
       width=4, height=12, dpi=300)

