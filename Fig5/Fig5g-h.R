
setwd('/data_hdd4/yuseong/cancer_aging/analysis/h.github/')

library(ggplot2)
library(ggpubr)
library(dplyr)
library(qs)
library(reshape2)
library(patchwork)


### load data
frac = qread('Fig2/data/2k-m//GSE193445_deconvolution_bisque_norm.qs')

frac = frac %>% as.data.frame() %>% tibble::rownames_to_column('Sample') %>% 
  mutate(Sample = gsub('Sample','Sample.', Sample))
head(frac)


meta = read.csv('Fig2/data/2k-m/GSE193445_META.csv')
head(meta)



#### Fig 5g - correlation between proportion of senescent cells
p_base <- ggplot(frac, aes(x = sen_malignant, y = sen_endothelial)) +
  geom_point(alpha = 0.8, size = 1.2, color = "black") +
  geom_smooth(method = "lm", color = "darkred", fill = "lightpink") +
  stat_cor(method = 'spearman', 
           label.x.npc = "left", label.y.npc = "top") +
  
  theme_bw() +
  labs(title = "GSE193445",
       x = "sen_malignant",
       y = "sen_endothelial")

p_base

ggsave('Fig5/figures/Fig5g_GSE193445_Bisque_sen_ratio_cor.pdf', width=6, height = 6, dpi=300)




#### Fig 5h - co-occurrence FRACTION permutation test 
set.seed(42)
N_PERM <- 10000

frac <- frac %>% as.data.frame() 

mal_prop  <- frac$sen_malignant
endo_prop <- frac$sen_endothelial

### Co-occurrence fraction = (Both / (Both + only-mal + only-endo))
cooccur_fraction <- function(mal, endo) {
  both      <- sum(mal & endo)
  only_mal  <- sum(mal & !endo)
  only_endo <- sum(!mal & endo)
  both / (both + only_mal + only_endo)
}

observed <- cooccur_fraction(mal_prop > 0, endo_prop > 0)
null_counts <- replicate(N_PERM, {
  mal_shuffled  <- sample(mal_prop) 
  endo_shuffled <- sample(endo_prop)
  cooccur_fraction(mal_shuffled > 0, endo_shuffled > 0)
})
p_val <- (1 + sum(null_counts >= observed)) / (N_PERM + 1)

## critical value: smallest achievable fraction f with P(null >= f) <= 0.05
tab <- sort(unique(null_counts))
crit <- min(tab[sapply(tab, function(f) mean(null_counts >= f) <= 0.05)])
crit_idx <- which(tab == crit)
crit_line <- if (crit_idx > 1) mean(c(tab[crit_idx - 1], crit)) else crit - 0.02

null_df <- as.data.frame(table(round(null_counts, 3)) / N_PERM)
colnames(null_df) <- c('n', 'freq')
null_df$n <- as.numeric(as.character(null_df$n))

sig_color <- if (p_val < 0.05) '#D62728' else '#1f9edb'
bar_width <- min(diff(sort(unique(null_df$n)))) * 0.85

p <- ggplot(null_df, aes(x = n, y = freq)) +
  geom_col(fill = 'grey75', width = bar_width) +
  geom_vline(xintercept = crit_line, linetype = 'dashed', color = 'black', linewidth = 0.6) +
  geom_vline(xintercept = observed, color = sig_color, linewidth = 0.9) +
  annotate('text', x = -Inf, y = Inf, label = paste0('observed = ', round(observed, 3)),
           hjust = -0.1, vjust = 1.8, size = 4.2, color = 'black') +
  annotate('text', x = -Inf, y = Inf, label = paste0('p = ', sprintf('%.4f', p_val)),
           hjust = -0.15, vjust = 3.4, size = 4.2, color = sig_color, fontface = 'bold') +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
  labs(title = 'Co-occurrence in GSE193445',
       x = 'co-occurrence fraction\n(Both / (Both + only-mal + only-endo))', y = 'probability frequency') +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(size = 15),
        plot.subtitle = element_text(size = 11),
        panel.grid.minor = element_blank(),
        legend.position = 'bottom')

legend_df <- data.frame(
  x = c(1, 2), y = c(1, 1),
  type = factor(c('observed samples (Sen malignant & Sen endothelial)', 'significance (p=0.05)'),
                levels = c('observed samples (Sen malignant & Sen endothelial)', 'significance (p=0.05)'))
)

p_legend <- ggplot(legend_df, aes(x, y, color = type, linetype = type)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = setNames(c(sig_color, 'black'), levels(legend_df$type)), name = NULL) +
  scale_linetype_manual(values = setNames(c('solid', 'dashed'), levels(legend_df$type)), name = NULL) +
  theme_void() +
  theme(legend.position = 'bottom', legend.text = element_text(size = 10)) +
  guides(color = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1))

final <- p / p_legend + plot_layout(heights = c(10, 1))
final

ggsave('Fig5/figures/Fig5h_bulkRT_sen_perm_cooccur_style_rowwise_both_bisque.pdf', final, width = 6, height = 5.5)



