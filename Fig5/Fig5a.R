
setwd('/data_hdd4/yuseong/cancer_aging/analysis/h.github/')

rm(list=ls())

library(CellChat)
library(patchwork)
library(qs)
library(dplyr)
library(Seurat)
library(tidyr)
library(stringr)
library(tibble)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)


### run CellChat

pat = qread('Fig2/data/2a-2d/patient_all_processed.qs')
pat$samples = pat$condition
table(pat$samples)

pat = subset(pat, subset = condition == 'Post')

cellchat <- createCellChat(object = pat, group.by = "celltype_f")

CellChatDB <- CellChatDB.human
head(CellChatDB.human$interaction)

# Users can now use this new database in CellChat analysis
load('Fig5/data/5a/CellChatDB.human_nichenetr_lrs.rda')

cellchat@DB <- db.new

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database

issue_genes <- c('C10orf99', 'CCL4L1', 'DEFB103B', 'DEFB106B', 'DEFB4A',
                 'DEFB4B', 'KLRA1P', 'LOC101059948', 'CGB8', 'LINC02210-CRHR1',
                 'RPSAP58', 'C19orf38')

filtered_interactions <- db.new$interaction %>%
  filter(!ligand %in% issue_genes, !receptor %in% issue_genes)

db.new$interaction <- filtered_interactions
cellchat@DB <- db.new

cellchat <- subsetData(cellchat)

future::plan("multisession", workers = 10) # do parallel
options(future.globals.maxSize = 20 * 1024 ^ 3)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat.oe = cellchat



ptm = Sys.time()

cellchat <- computeCommunProb(cellchat, type = "triMean", trim=0.1)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

qsave(cellchat, 'Fig5/data/5a/CellChat_results2_post_aggregated_com.qs', preset='fast')


cellchat = qread('Fig5/data/5a/CellChat_results2_post_aggregated_com.qs')


### Fig 5a - sen mal interactions
targets = c("MuscleCell","PlasmaCell","Bcell","Myeloid","T.NKcell","non_malignant",
            "Endothelial","Fibroblast","malignant" )
n_tar = length(targets)

all_sectors <- c("sen_malignant", targets)
n_sectors <- length(all_sectors)

circos.clear()
circos.par(gap.after = rep(6, n_sectors), start.degree = 35)

mock_values <- cellchat@net$weight['sen_malignant', targets]
color_palette <- colorRampPalette(c("#2166ac", "#FFFFBF", "#b2182b"))(100)

sector_colors <- c("sen_malignant" = "#333333")

min_val <- min(mock_values)
max_val <- max(mock_values)

for(cell_name in names(mock_values)) {
  current_val <- mock_values[cell_name]

  if (max_val == min_val) {
    color_idx <- 50
  } else {
    color_idx <- round((current_val - min_val) / (max_val - min_val) * 99) + 1
  }
  sector_colors[cell_name] <- color_palette[color_idx]
}


pdf("Fig5/figures/sen_mal_CellChat_circleplt.pdf", width = 6, height = 6)

circos.clear()
circos.par(gap.after = rep(6, n_sectors), start.degree = 35)

sector_widths <- c("sen_malignant" = 8, setNames(rep(4, n_tar), targets))
initialize_df <- data.frame(factor(all_sectors, levels = all_sectors), x = 0, y = 1)
circos.initialize(factors = initialize_df[,1], xlim = c(0, 1), sector.width = sector_widths)

circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.name = get.cell.meta.data("sector.index")
  if(sector.name == "sen_malignant") {
    circos.points(x = 0.5, y = 0.5, pch = 16, col = sector_colors[sector.name], cex = 5.5)
    circos.text(x = 0.5, y = 0.5, labels = "Sen_Mal", col = "white", cex = 0.7, font = 2)
  } else {
    current_weight <- mock_values[sector.name]
    calculated_cex <- (current_weight / max_val) * 5.0
    if(calculated_cex < 0.7) { calculated_cex <- 0.7 }

    circos.points(x = 0.5, y = 0.5, pch = 16, col = sector_colors[sector.name], cex = calculated_cex)
    circos.points(x = 0.5, y = 0.5, pch = 1, col = "black", cex = calculated_cex)

    text_y_margin <- 1.4 + (calculated_cex * 0.15)
    circos.text(x = 0.5, y = text_y_margin, labels = sector.name,
                niceFacing = TRUE, adj = c(0, 0.5), cex = 0.75, font = 2)
  }
}, bg.border = NA, track.height = 0.1)

for(i in 1:length(targets)) {
  cell <- targets[i]
  current_weight <- mock_values[cell]
  calculated_w <- (current_weight / max_val) * 10
  if(calculated_w < 0.8) { calculated_w <- 0.8 }

  link_color <- adjustcolor("#b5b5b5", alpha.f = 0.45)

  circos.link(
    sector.index1 = "sen_malignant", point1 = 0.5,
    sector.index2 = cell, point2 = 0.5,
    col = link_color, border = NA, lwd = calculated_w
  )
}


color_fun = colorRamp2(seq(min_val, max_val, length.out = 4),
                       colorRampPalette(c("#2166ac", "#FFFFBF", "#b2182b"))(4))

lgd_color = Legend(
  title = "interaction weight",
  col_fun = color_fun,
  title_gp = gpar(fontsize = 10, fontface = "bold"),
  labels_gp = gpar(fontsize = 8, fontface = "bold"),
  grid_width = unit(5, "mm"),
  legend_height = unit(30, "mm")
)


legend_breaks <- c(1, 6, 12, 18)
legend_sizes  <- (legend_breaks / max_val) * 5.0
legend_sizes[legend_sizes < 0.7] <- 0.7

lgd_size = Legend(
  title = "Node Size (Weight)",
  type = "points",
  pch = 16,
  size = unit(legend_sizes * 2.5, "mm"),
  labels = as.character(legend_breaks),

  title_gp = gpar(fontsize = 10, fontface = "bold"),
  labels_gp = gpar(fontsize = 8, fontface = "bold"),
  grid_height = unit(8, "mm")
)

pack_legends = packLegend(lgd_color, lgd_size, direction = "vertical")

draw(pack_legends, x = unit(0.88, "npc"), y = unit(0.5, "npc"), just = c("center", "center"))

circos.clear()
dev.off()


### Fig 5a - sen endo interactions
targets2 = c("MuscleCell","PlasmaCell","Bcell","Myeloid","T.NKcell","non_malignant",
             "Endothelial","Fibroblast","malignant" )
n_tar2 = length(targets2)

all_sectors2 <- c("sen_endothelial", targets2)
n_sectors2 <- length(all_sectors2)

circos.clear()
circos.par(gap.after = rep(6, n_sectors2), start.degree = 35)

mock_values2 <- cellchat@net$weight['sen_endothelial', targets2]
color_palette <- colorRampPalette(c("#2166ac", "#FFFFBF", "#b2182b"))(100)

sector_colors2 <- c("sen_endothelial" = "#333333")

min_val2 <- min(mock_values2)
max_val2 <- max(mock_values2)

for(cell_name in names(mock_values2)) {
  current_val <- mock_values2[cell_name]

  if (max_val2 == min_val2) {
    color_idx <- 50
  } else {
    color_idx <- round((current_val - min_val2) / (max_val2 - min_val2) * 99) + 1
  }
  sector_colors2[cell_name] <- color_palette[color_idx]
}

pdf("Fig5/figures/Fig5a_sen_endo_CellChat_circleplt.pdf", width = 6, height = 6)

circos.clear()
circos.par(gap.after = rep(6, n_sectors2), start.degree = 35)

sector_widths2 <- c("sen_endothelial" = 8, setNames(rep(4, n_tar2), targets2))
initialize_df2 <- data.frame(factor(all_sectors2, levels = all_sectors2), x = 0, y = 1)
circos.initialize(factors = initialize_df2[,1], xlim = c(0, 1), sector.width = sector_widths2)

circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.name = get.cell.meta.data("sector.index")
  if(sector.name == "sen_endothelial") {
    circos.points(x = 0.5, y = 0.5, pch = 16, col = sector_colors2[sector.name], cex = 5.5)
    circos.text(x = 0.5, y = 0.5, labels = "Sen_Endo", col = "white", cex = 0.7, font = 2)
  } else {
    current_weight <- mock_values2[sector.name]
    calculated_cex <- (current_weight / max_val2) * 5.0
    if(calculated_cex < 0.7) { calculated_cex <- 0.7 }

    circos.points(x = 0.5, y = 0.5, pch = 16, col = sector_colors2[sector.name], cex = calculated_cex)
    circos.points(x = 0.5, y = 0.5, pch = 1, col = "black", cex = calculated_cex)

    text_y_margin <- 1.4 + (calculated_cex * 0.15)
    circos.text(x = 0.5, y = text_y_margin, labels = sector.name,
                niceFacing = TRUE, adj = c(0, 0.5), cex = 0.75, font = 2)
  }
}, bg.border = NA, track.height = 0.1)

for(i in 1:length(targets2)) {
  cell <- targets2[i]
  current_weight <- mock_values2[cell]
  calculated_w <- (current_weight / max_val2) * 10
  if(calculated_w < 0.8) { calculated_w <- 0.8 }
  
  link_color <- adjustcolor("#b5b5b5", alpha.f = 0.45)

  circos.link(
    sector.index1 = "sen_endothelial", point1 = 0.5,
    sector.index2 = cell, point2 = 0.5,
    col = link_color, border = NA, lwd = calculated_w
  )
}


color_fun2 = colorRamp2(seq(min_val2, max_val2, length.out = 4),
                        colorRampPalette(c("#2166ac", "#FFFFBF", "#b2182b"))(4))

lgd_color = Legend(
  title = "interaction weight",
  col_fun = color_fun2,
  title_gp = gpar(fontsize = 10, fontface = "bold"),
  labels_gp = gpar(fontsize = 8, fontface = "bold"),
  grid_width = unit(5, "mm"),
  legend_height = unit(30, "mm")
)

legend_breaks <- c(1, 6, 12, 18)
legend_sizes  <- (legend_breaks / max_val2) * 5.0
legend_sizes[legend_sizes < 0.7] <- 0.7

lgd_size = Legend(
  title = "Node Size (Weight)",
  type = "points",
  pch = 16,
  size = unit(legend_sizes * 2.5, "mm"),
  labels = as.character(legend_breaks),

  title_gp = gpar(fontsize = 10, fontface = "bold"),
  labels_gp = gpar(fontsize = 8, fontface = "bold"),
  grid_height = unit(8, "mm")
)

pack_legends = packLegend(lgd_color, lgd_size, direction = "vertical")

draw(pack_legends, x = unit(0.88, "npc"), y = unit(0.5, "npc"), just = c("center", "center"))

circos.clear()
dev.off()



