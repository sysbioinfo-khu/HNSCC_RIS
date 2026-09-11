
setwd('/data_hdd4/yuseong/cancer_aging/analysis/h.github/Fig2')

rm(list=ls())

library(qs)     
library(Seurat)          
library(omnideconv) 
library(dplyr)             
library(tibble)         
library(tidyr)        
library(reshape2)    
library(ggplot2)         
library(ggpubr)  
library(survival)         
library(survminer)  


# SET CORES TO USE
library(BiocParallel)
library(RhpcBLASctl)
blas_set_num_threads(35)
omp_set_num_threads(35)



## cell type deconvolution via Bisque
# load sc for reference
pat = qread('data/2a-2d/patient_seurat_annotated_FINAL+copykat+senanno.qs')
sc_counts <- as.matrix(Seurat::GetAssayData(pat, assay = "RNA", layer = "counts"))
cell_types <- as.character(pat$celltype_f)
batch_ids <- as.character(pat$Sample)

# load bulk count data
bulk_counts = read.csv('data/2k-m/GSE193445_counts.csv', row.names=1)

# run deconvolution
frac <- omnideconv::deconvolute(bulk_gene_expression = bulk_counts,
                                         method = 'bisque',
                                         single_cell_object = sc_counts, 
                                         batch_ids = batch_ids,
                                         cell_type_annotations = cell_types, 
                                         normalize_results = TRUE
)

qsave(frac, 'data/2k-m/GSE193445_deconvolution_bisque_norm.qs', preset='fast')


### Fig 2k - comparison of senescence phenotypes between pre- and post-chemoradiotherapy samples
meta = read.csv('data/2k-m/GSE193445_META.csv')
meta = meta %>% dplyr::rename('group' = Response_final)
head(meta)

## Fig 2k left panel
frac = qread('data/2k-m/GSE193445_deconvolution_bisque_norm.qs')

frac = frac %>% as.data.frame() %>% rownames_to_column('Sample') %>% 
  mutate(Sample = gsub('Sample','Sample.', Sample))
frac$sen_total = frac$sen_endothelial + frac$sen_malignant
head(frac)

frac_long <- melt(frac)
colnames(frac_long) <- c("Sample", "celltype", "Score")

frac_long = left_join(frac_long, meta)
frac_long$condition = factor(frac_long$condition,  levels = c("pre", "post"))                        

plot_df <- frac_long %>%
  filter(celltype == 'sen_total') %>%
  group_by(celltype, condition, group) %>%
  summarise(
    mean  = mean(Score),
    sem   = sd(Score) / sqrt(n()),
    n     = n(),
    .groups = "drop"
  )

pval_df <- frac_long %>%
  filter(celltype == 'sen_total', condition == "post") %>%
  group_by(celltype) %>%
  summarise(
    p = t.test(Score ~ group)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    label = case_when(
      p <= 0.0001 ~ '****',
      p <= 0.001 ~ "***",
      p <= 0.01  ~ "**",
      p <= 0.05  ~ "*",
      TRUE       ~ "ns"
    ),
    condition = "post",
    x = 2,
    y = 0.05        
  )

ggplot(plot_df, aes(x = condition, y = mean, group = group, color = group)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.1) +
  facet_wrap(~ celltype, scales = "free_y", ncol = 3) +
  geom_text(data = pval_df,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 5) +
  scale_color_manual(values = c("CR" = "#0099FF", "NR" = "#d03b3b")) +
  labs(y = "Scaled GSVA score", x = NULL, color = "Response") +
  ylim(0, 0.05) +
  theme_bw()

ggsave('figures/Fig2k_GSE193445_slope_sen_bisque.pdf', width=4, height =5)


## Fig 2k other panels
gsva.es = read.csv('data/2k-m/GSE193445_GSVA_pathwayscores.csv', row.names=1)

gsva_scaled <- t(apply(t(gsva.es), 1, function(x) (x - min(x)) / (max(x) - min(x))))
gsva_scaled = as.data.frame(t(as.data.frame(gsva_scaled))) %>% rownames_to_column('Sample') %>% filter(Sample %in% sel.p)

gsva_long <- melt(gsva_scaled)
colnames(gsva_long) <- c("Sample", "cluster", "Score")

gsva_long = left_join(gsva_long, meta)

gsva_long$condition = factor(gsva_long$condition,  levels = c("pre", "post"))                        

ggplot(gsva_long[gsva_long$cluster %in% names(gn_ls),], aes(x = condition, y = Score, group = group, color = group)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  stat_summary(fun = mean, geom = "point", size = 2.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
  facet_wrap(~ cluster, scales = "free_y", ncol = 4) +
  stat_compare_means(
    data = subset(gsva_long, condition == "post" & cluster %in% names(gn_ls)),
    aes(x = condition, y = Score, group = group),
    method = "t.test",   
    label = "p.signif"          
  )+
  scale_color_manual(values = c("CR" = "#0099FF", "NR" = "#d03b3b")) +
  labs(y = "Scaled GSVA score", x = NULL, color = "Response") +
  ylim(0,1) +
  theme_bw()


ggsave('figures/Fig2k_GSE193445_slope_GSVA.pdf', width=10, height =5)




### Fig 2l - correlation between proportion of RIS popoulation and pathway scores
gsva_frac = inner_join(gsva_scaled, frac)
gf.merge = inner_join(gsva_frac, meta)

rt.dat = gf.merge[gf.merge$condition == 'post',]
length(unique(rt.dat$Sample))

plot_df <- rt.dat %>%
  select(Sample, sen_total, SEN_MAYO, EMT, CELLULAR.RESPONSE.TO.UV, M3, group) %>%
  pivot_longer(cols = c(SEN_MAYO,EMT, CELLULAR.RESPONSE.TO.UV, M3),
               names_to = "phenotype", values_to = "score")
head(plot_df)

ggplot(plot_df, aes(x = sen_total, y = score)) +
  geom_point(aes(color = group), size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = F, color = "black", linewidth = 0.5) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top", size = 3.5) +
  facet_wrap(~ phenotype, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("CR" = "#0099FF", "NR" = "#d03b3b")) +
  ylim(0, 1.3) +
  labs(x = "proportion of senescence cells", y = "score") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave('figures/Fig2l_GSE193445_cor_bisque.pdf', width=8, height =5)




### Fig 2m - Kaplan-Meier analysis in the RIS population
gsva_frac = inner_join(gsva_scaled, frac)
gsva_frac_long <- melt(gsva_frac)
colnames(gsva_frac_long) <- c("Sample", "cluster","Score")

gf.merge_long = inner_join(gsva_frac_long, meta)

rt.dat_post = gf.merge_long[gf.merge_long$cluster == 'sen_total' & gf.merge_long$condition == 'post',]
rt.dat_post$sample = gsub('\\.', '', rt.dat_post$Sample)

cof <- median(rt.dat_post$Score)
cof

rt.dat_post$sen_g <- ifelse(rt.dat_post['Score'] > cof, 'high', 'low')
table(rt.dat_post$sen_g)

sm.high = rt.dat_post[rt.dat_post$sen_g == 'high',]$sample

surv.rt.dat = rt.dat_post %>%
  dplyr::select(sample, Recurrence, DFS, sen_g) %>%
  mutate(Recurrence = case_when(Recurrence == 'yes'~ 1,
                                Recurrence == 'no' ~ 0)) %>% 
  distinct()

surv_test <- survdiff(Surv(DFS, Recurrence) ~ sen_g, data = surv.rt.dat)
pval <- 1 - pchisq(surv_test$chisq, df = length(surv_test$n) - 1)
fit <- survfit(Surv(DFS, Recurrence) ~ sen_g, data=surv.rt.dat)

gp_cnts <- table(surv.rt.dat$sen_g)
gp_labels <- paste0(names(gp_cnts), " (n=", gp_cnts, ")")

p2 <- ggsurvplot(fit, data = surv.rt.dat,
                 pval=T, xlab = 'DFS (months)',
                 ylab = 'Recurrence-free probability',
                 legend.labs = gp_labels,
                 # risk.table = TRUE,
                 # break.time.by = 12,
                 title = 'GSE193445',
                 font.main = c(20),
                 font.legend = 15,
                 palette = c("darkorange2", "lightblue4"))

p2

ggsave('figures/Fig2m_GSE193445_sen_surv_recur.pdf', width=8, height=6)
