
rm(list=ls())

setwd('/data_hdd4/yuseong/cancer_aging/analysis/h.github/Fig3')


library(jsonlite)
library(data.table)
library(igraph)


tree_paths <- read_json('Fig3/data/3e-f/alpaca/input/HyPR.HN.02/tree_paths.json')
tree_edges = tree_paths

### PLOT TREE OF PRE
de_emphasize_clones <- c('C6')
common_clones = c('C1','C3','C7')

extract_edges <- function(path) {
  if (length(path) < 2) return(NULL)
  cbind(path[-length(path)], path[-1])
}

all_edges <- do.call(rbind, lapply(tree_paths, extract_edges)) %>%
  as.data.frame() %>%
  distinct()

colnames(all_edges) <- c("from", "to")

all_edges$from = gsub('clone','C', all_edges$from)
all_edges$to = gsub('clone','C', all_edges$to)

g <- graph_from_data_frame(all_edges, directed = TRUE)

v_names <- V(g)$name

trunk = intersect(tree_edges[[1]], tree_edges[[2]])
trunk = unlist(intersect(trunk,tree_edges[[3]]))
trunk = gsub('clone','C',trunk)
l <- layout_as_tree(g, root = trunk)

final_vcol <- ifelse(v_names %in% de_emphasize_clones, "white", 
                     ifelse(v_names %in% common_clones, '#AEEEEE', 
                            ifelse(v_names == 'C4', '#A2CD5A', 
                                   ifelse(v_names == 'C8', '#27408B', '#1E90FF'))))

final_vsize <- ifelse(v_names %in% de_emphasize_clones, 7, 15)
final_vframe <- ifelse(v_names %in% de_emphasize_clones, "grey80", "black")
final_vlwd <- ifelse(v_names %in% de_emphasize_clones, 1, 20)

edge_list <- as_edgelist(g)
ecol <- apply(edge_list, 1, function(x) {
  if (any(x %in% de_emphasize_clones)) "grey80" else "black"
})

ewid <- apply(edge_list, 1, function(x) {
  if (any(x %in% de_emphasize_clones)) 1 else 2.5
})

pdf("figures/Fig3e_WES_treeplot_pre.pdf", width = 5, height = 5)
plot.igraph(g
            , layout = l
            , vertex.color = final_vcol
            , vertex.size = final_vsize
            , vertex.frame.color = final_vframe
            , vertex.lwd = final_vlwd
            , edge.color = ecol        
            , edge.width = ewid        
            , vertex.label = c(rep('', 8))
            , vertex.label.color = "black"
            , vertex.label.family = "sans"
            , vertex.label.font = 2
            , edge.arrow.size = 0 
            , main = '')

dev.off()




### PLOT TREE OF POST
tree_edges = tree_paths

de_emphasize_clones <- c('C2','C8', 'C9')
common_clones = c('C1','C3','C7')

v_names <- V(g)$name

final_vcol <- ifelse(v_names %in% de_emphasize_clones, "white", 
                     ifelse(v_names %in% common_clones, '#AEEEEE', 
                            ifelse(v_names == 'C4', '#A2CD5A', '#FF6347')))
final_vsize <- ifelse(v_names %in% de_emphasize_clones, 7, 15)
final_vframe <- ifelse(v_names %in% de_emphasize_clones, "grey80", "black")
final_vlwd <- ifelse(v_names %in% de_emphasize_clones, 1, 20)

edge_list <- as_edgelist(g)
ecol <- apply(edge_list, 1, function(x) {
  if (any(x %in% de_emphasize_clones)) "grey80" else "black"
})

ewid <- apply(edge_list, 1, function(x) {
  if (any(x %in% de_emphasize_clones)) 1 else 2.5
})


pdf("figures/Fig3e_WES_treeplot_post.pdf", width = 5, height = 5)
plot.igraph(g
            , layout = l
            , vertex.color = final_vcol
            , vertex.size = final_vsize
            , vertex.frame.color = final_vframe
            , vertex.lwd = final_vlwd
            , edge.color = ecol   
            , edge.width = ewid       
            , vertex.label = c(rep('', 8))
            , vertex.label.color = "black"
            , vertex.label.family = "sans"
            , vertex.label.font = 2
            , edge.arrow.size = 0   
            , main = '')

dev.off()


