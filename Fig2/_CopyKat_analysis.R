
setwd('/data_hdd4/yuseong/cancer_aging/analysis/h.github/Fig2')

library(Seurat)
library(SeuratObject)
library(qs)
library(copykat)


#### load data
pat = qread('data/2a-2d/patient_all_processed.qs')
colnames(pat) <- gsub('-', '.', colnames(pat))
pat$condition = factor(pat$condition, levels = c('Pre', 'On', 'Post'))

cnts = GetAssayData(pat, slot="counts")

sample_ids <- sapply(strsplit(colnames(cnts), "\\."), function(x) {paste(x[1:3], collapse='.')})
cnts_list <- split(colnames(cnts), sample_ids)
cnts_by_sample <- lapply(cnts_list, function(cells) cnts[, cells, drop=FALSE])


### RUN CopyKAT 
## Liu, Z.L., Meng, X.Y., Bao, R.J. et al. 
## Single cell deciphering of progression trajectories of the tumor ecosystem in head and neck cancer.
## Nat Commun 15, 2595 (2024).

# id.type=“S”, ngene.chr=1, win.size=25, KS.cut=0.1, sam.name=“test”, distance=“euclidean”,
# n.cores=20,output.seg=“FLASE”.

for (i in 1:4) {
  sam = names(cnts_by_sample)[[i]]
  outp = paste0('data/others/copykat/',sam)
  dir.create(outp)
  
  setwd(outp)
  
  cnts = cnts_by_sample[[i]]
  nonzero_gene_idx <- Matrix::rowSums(cnts) > 0
  cnts_filtered <- cnts[nonzero_gene_idx, ]
  
  copykat.test <- copykat(rawmat=cnts_filtered, 
                          id.type="S", 
                          ngene.chr=1, 
                          win.size=25, 
                          KS.cut=0.1, 
                          sam.name=sam, 
                          distance="euclidean", 
                          output.seg="FLASE", 
                          plot.genes="TRUE", 
                          genome="hg20",
                          n.cores=20)
  
  print(paste0(sam, ' is finished!'))
  print('++++++++++++++++++++++++++++++++++++++')
}


### OUTPUT
sam_d = list.dirs(path='data/others/copykat', full.names = T)
sam_d = grep('/Hy', sam_d, value = T)

pred = data.frame()
for (sam in sam_d) {
  pred1 = read.table(paste0(sam, '/test_copykat_prediction.txt'), sep='\t', header=T)
  pred = rbind(pred, pred1)
}

ck = pred$copykat.pred
names(ck) = pred$cell.names

pat$copykat_nonref = ck

epi.so = subset(pat, subset = CellType == 'Epithelial')

table(epi.so$copykat_nonref)




