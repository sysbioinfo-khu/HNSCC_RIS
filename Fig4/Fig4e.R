setwd('/data_hdd4/yuseong/cancer_aging/analysis/h.github/Fig4')

rm(list=ls())

library(readxl)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)


tf_list <- list(
  POU2F2 = "POU2F2_M07972",
  EGR1 = "EGR1_M07927",
  FLI1 = "FLI1_M08111",
  NFATC2 = "NFATC2_M01511"
)

all_tf_data <- lapply(names(tf_list), function(tf_id) {
  file_path <- paste0('data/4e/tobias/', tf_list[[tf_id]], '/', tf_list[[tf_id]], '_overview.xlsx')
  
  df <- read_excel(file_path) %>%
    filter(!(`0GY_score` == 0 & `8GY_score` == 0)) %>%
    filter(abs(`0GY_8GY_log2fc`) > 2) %>%
    dplyr::select(peak_id, `0GY_score`, `8GY_score`) %>%
    group_by(peak_id) %>%
    summarise(WT = mean(`0GY_score`, na.rm = TRUE),
              RT = mean(`8GY_score`, na.rm = TRUE)) %>%
    ungroup()
  
  df$TF <- tf_id
  return(df)
}) %>% bind_rows()


all_tf_long <- melt(all_tf_data, id.vars = c("peak_id", "TF"), 
                    variable.name = "Condition", value.name = "Binding_score")

all_tf_long$TF = factor(all_tf_long$TF, levels = c('POU2F2','EGR1','FLI1','NFATC2'))


stats_df <- all_tf_long %>%
  arrange(TF, peak_id, Condition) %>%
  group_by(TF) %>%
  summarise(
    p_value = t.test(
      Binding_score[Condition == "WT"], 
      Binding_score[Condition == "RT"], 
      paired = TRUE
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(p.adj = p.adjust(p_value, method = "BH")) %>%
  mutate(p.signif = case_when(
    p.adj < 0.001 ~ "***",
    p.adj < 0.01  ~ "**",
    p.adj < 0.05  ~ "*",
    TRUE          ~ "ns"
  ))

stats_df


ggplot(all_tf_long, aes(x = TF, y = Binding_score, fill = Condition)) +
  stat_summary(fun = "mean", geom = "bar", 
               position = position_dodge(width = 0.8), 
               width = 0.7, color = "black") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", 
               position = position_dodge(width = 0.8), 
               width = 0.2) +
  geom_text(data = stats_df, aes(x = TF, y = 0.6, label = p.signif), 
            inherit.aes = FALSE, vjust = -0.2, fontface = "bold") +
  scale_fill_manual(values = c("WT" = "#999999", "RT" = "#E64B35FF")) + 
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text  = element_text(color = "black"),
    legend.title = element_text(face = "bold")
  ) +
  labs(x = "Transcription Factors", y = "Binding Score")


ggsave('figures/Fig4e_ATAC_TFbindingscore_barplt.pdf',
       width=8,height = 6, dpi=300)




