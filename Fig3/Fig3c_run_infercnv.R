setwd('/data_hdd4/yuseong/cancer_aging/analysis/h.github/')

library(Seurat)
library(SeuratObject)
library(dplyr)
library(stringr)
library(Matrix)
library(infercnv)
library(qs)


#############################
###  MAKE INFERCNV INPUTS ###
#############################
infercnv_run = function(pat, pat.sc) {
  mal = subset(pat.sc, subset = CellType == 'Epithelial')
  mal2 = subset(mal, subset = grep(pat, Sample))
  unique(mal2$Sample)
  
  outdir = paste0('Fig3/data/3c/infercnv/', pat)
  dir.create(outdir)
  
  ## INPUT1 : SAMPLE ANNOTATION
  anno = as.data.frame(list('cnv' = mal2$copykat_nonref))
  head(anno)
  
  anno$celltype = ifelse(anno$cnv == 'aneuploid', 'Malignant','Normal')
  head(anno)
  
  anno$Sample = sapply(strsplit(rownames(anno), '_'), function(x) {x[1]})
  head(anno)
  unique(anno$Sample)
  
  anno$cell = gsub('-','.', rownames(anno))
  anno$celltype = ifelse(anno$celltype=='Malignant', paste0(anno$celltype,'_', anno$Sample), anno$celltype)
  rownames(anno) = NULL
  anno_f = anno[c('cell', 'celltype')]
  
  head(anno_f)
  
  unique(anno$celltype)
  
  write.table(anno_f, file.path(outdir, "pat.epi_ck.sample.annot.txt"), sep='\t', row.names=F, col.names=F, quote=F)
  
  ## INPUT2 : COUNT MATRIX
  counts_matrix = GetAssayData(mal2, slot="counts")
  colnames(counts_matrix) = gsub('-','.', colnames(counts_matrix))
  counts_matrix[1:5, 1:5]
  
  gene_sums <- rowSums(counts_matrix)
  counts_matrix <- counts_matrix[gene_sums > 0, ]
  
  ## INPUT3 : GENOMIC POSITION ANNOTATION
  ## make genomic position file
  # python /data_hdd4/yuseong/cancer_aging/analysis/h.github/Fig3/analysis/gtf_to_position_file.py \
  # --attribute_name gene_name \
  # /data_hdd4/yuseong/reference/hg38/gencode.v48.primary_assembly.annotation.gtf \
  # /data_hdd4/yuseong/cancer_aging/analysis/h.github/Fig3/data/3c/infercnv/gencode.v48.primary_assembly.annotation.gtf.gen_pos.txt
  
  
  ######################
  ###  RUN INFERCNV  ###
  ######################
  # create the infercnv object
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                      annotations_file=file.path(outdir, "pat.epi_ck.sample.annot.txt"),
                                      delim="\t",
                                      gene_order_file="Fig3/data/3c/infercnv/gencode.v48.primary_assembly.annotation.gtf.gen_pos.txt",
                                      ref_group_names=c("Normal"))
  
  str(infercnv_obj)
  
  # perform infercnv operations to reveal cnv signal
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                               out_dir=outdir,  # dir is auto-created for storing outputs
                               denoise=T,
                               HMM=T,
                               num_threads = 20,
                               cluster_by_groups = FALSE
  )
}

pat.sc = qread('Fig2/data/2a-2d/patient_all_processed.qs')
pat.sc$SAMPLE = sapply(strsplit(pat.sc$Sample, split ='-'), function(x) {paste0(x[1:3], collapse='-')})
unique(pat.sc$SAMPLE)

epi = subset(pat.sc, subset = CellType == 'Epithelial')
qsave(epi, 'Fig3/data/3c/patient_epithelial_subset.qs', preset='fast')


for (sam in unique(pat.sc$SAMPLE)) {
  infercnv_run(sam, pat.sc)
  print(sam)
}




