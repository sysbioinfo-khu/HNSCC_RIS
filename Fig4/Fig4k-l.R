setwd('/data_hdd4/yuseong/cancer_aging/analysis/h.github/Fig4')

rm(list=ls())

library(ggplot2)
library(patchwork)
library(dplyr)
library(readxl)


### plot functions

fmt_p <- function(p) {
  if (p < 0.001) "P < 0.001" else sprintf("P = %.3f", p)
}

make_4cond_plot <- function(summary_df, points_df, p_left, p_right, y_lab) {
  
  ytop <- max(summary_df$mean + summary_df$se, points_df$value) * 1.35
  y1 <- max(c(summary_df$mean[summary_df$pos %in% c(1, 2)] +
                summary_df$se[summary_df$pos %in% c(1, 2)],
              points_df$value[points_df$pos %in% c(1, 2)])) * 1.30
  y2 <- max(c(summary_df$mean[summary_df$pos %in% c(3, 4)] +
                summary_df$se[summary_df$pos %in% c(3, 4)],
              points_df$value[points_df$pos %in% c(3, 4)])) * 1.20
  
  main <- ggplot(summary_df, aes(x = pos, y = mean, fill = siRNA)) +
    geom_col(width = 0.6, color = "black") +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.15) +
    geom_point(data = points_df, aes(x = pos, y = value), shape = 21,
               fill = "white", color = "black", size = 2, inherit.aes = FALSE) +
    annotate("segment", x = 1, xend = 2, y = y1, yend = y1, linewidth = 0.4) +
    annotate("text", x = 1.5, y = y1 * 1.10, label = fmt_p(p_left), size = 3.5) +
    annotate("segment", x = 3, xend = 4, y = y2, yend = y2, linewidth = 0.4) +
    annotate("text", x = 3.5, y = y2 * 1.08, label = fmt_p(p_right), size = 3.5) +
    scale_fill_manual(values = c("Ctrl" = "#9ab1bb", "KD" = "#e36c6e"), guide = "none") +
    scale_x_continuous(limits = c(0.4, 4.6), breaks = NULL, expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)), limits = c(0, ytop)) +
    labs(x = NULL, y = y_lab) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text = element_text(color = "black"),
      plot.margin = margin(t = 10, r = 10, b = 2, l = 5)
    )
  
  table_df <- data.frame(
    pos = rep(1:4, 2),
    row = factor(rep(c("IR (7.5 Gy)", "POU2F2 siRNA"), each = 4),
                 levels = c("IR (7.5 Gy)", "POU2F2 siRNA")),
    label = c("-", "-", "+", "+",
              "-", "+", "-", "+")
  )
  
  tbl <- ggplot(table_df, aes(x = pos, y = row, label = label)) +
    geom_text(size = 3.8) +
    scale_x_continuous(limits = c(0.4, 4.6), expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_text(color = "black", hjust = 1),
      panel.grid = element_blank(),
      plot.margin = margin(t = 0, r = 10, b = 5, l = 5)
    )
  
  main / tbl + plot_layout(heights = c(5, 1.3))
}



### Fig 4k - SA-B-gal staining results in ctrl siRNA and POU2F2 siRNA
raw <- read_excel("data/4k-l/POU2F2siRNA_SA_beta_gal_raw_count.xlsx")

raw$Experiment <- factor(raw$Experiment)
raw$Group <- factor(raw$Group, levels = c("Ctrl siRNA", "POU2F2 siRNA"))
raw$Treatment <- factor(raw$Treatment, levels = c("NR", "IR"))

# Pooled ratio per biological replicate: Sigma(Positive_cells) / Sigma(Total_cells)
df_bio <- raw %>%
  group_by(Experiment, Group, Treatment) %>%
  summarise(Bio_mean = sum(Positive_cells) / sum(Total_cells) * 100, .groups = "drop") %>%
  mutate(siRNA = ifelse(Group == "Ctrl siRNA", "Ctrl", "KD"),
         pos = case_when(
           Treatment == "NR" & siRNA == "Ctrl" ~ 1,
           Treatment == "NR" & siRNA == "KD"   ~ 2,
           Treatment == "IR" & siRNA == "Ctrl" ~ 3,
           Treatment == "IR" & siRNA == "KD"   ~ 4
         ))

summary_df <- df_bio %>%
  group_by(pos, siRNA) %>%
  summarise(mean = mean(Bio_mean), se = sd(Bio_mean) / sqrt(n()), .groups = "drop")

points_df <- df_bio %>% transmute(pos, value = Bio_mean, siRNA)

levels(df_bio$Group) # "Ctrl siRNA"   "POU2F2 siRNA"

p_left  <- t.test(Bio_mean ~ Group, data = filter(df_bio, Treatment == "NR"),  alternative = "greater")$p.value
p_right <- t.test(Bio_mean ~ Group, data = filter(df_bio, Treatment == "IR"), alternative = "greater")$p.value

p <- make_4cond_plot(summary_df, points_df, p_left, p_right,
                     y_lab = expression(paste("SA-", beta, "-gal positive cells (%)")))
p

ggsave("figures/Fig4k_siRNA_SA_beta_gal_4cond.pdf", plot = p, width = 4, height = 6, units = "in", dpi = 300)


### Fig 4l - sorted senescent cells in ctrl siRNA and POU2F2 siRNA
facs.df = read.csv("data/4k-l/POU2F2siRNA_sorted_senescent_cells_result.csv")

summary_df <- facs.df %>%
  group_by(pos, siRNA) %>%
  summarise(mean = mean(fraction), se = sd(fraction) / sqrt(n()), .groups = "drop")

points_df <- facs.df %>% transmute(pos, value = fraction, siRNA)

p_left <- t.test(
  facs.df$fraction[facs.df$radiation == "NR" & facs.df$siRNA == "Ctrl"],
  facs.df$fraction[facs.df$radiation == "NR" & facs.df$siRNA == "KD"],
  paired = TRUE, alternative = "greater"
)$p.value
p_right <- t.test(
  facs.df$fraction[facs.df$radiation == "IR" & facs.df$siRNA == "Ctrl"],
  facs.df$fraction[facs.df$radiation == "IR" & facs.df$siRNA == "KD"],
  paired = TRUE, alternative = "greater"
)$p.value

p2 <- make_4cond_plot(summary_df, points_df, p_left, p_right,
                     y_lab = "Senescent cell fraction (%)")
p2

ggsave("figures/Fig4l_siRNA_Flow_senescence_4cond.pdf", plot = p2,  width = 4, height = 6, units = "in", dpi = 300)

