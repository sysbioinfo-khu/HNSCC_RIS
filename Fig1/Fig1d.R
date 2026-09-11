
setwd('/data_hdd4/yuseong/cancer_aging/analysis/h.github/Fig1')

rm(list=ls())

library(ggplot2)
library(dplyr)
library(patchwork)  

### Fig 1d - SA-Beta-gal staining results

df_3g = read.csv('data/1d/EXP_RTcellline_SA-B-Gal_prop.csv')

summary_3g <- df_3g %>%
  group_by(Group) %>%
  summarise(
    Mean = mean(Ratio),
    SEM  = sd(Ratio) / sqrt(n())
  )

make_donut <- function(group_name, color_pos) {
  stats <- summary_3g %>% filter(Group == group_name)
  mean_val <- round(stats$Mean, 1)
  sem_val  <- round(stats$SEM, 1)
  
  pie_df <- data.frame(
    Category = factor(c("Positive", "Negative"),
                      levels = c("Positive", "Negative")),
    Value = c(mean_val, 100 - mean_val)
  )
  
  ggplot(pie_df, aes(x = 2, y = Value, fill = Category)) +
    geom_bar(stat = "identity", width = 1,
             color = "white", linewidth = 0.8) +
    coord_polar(theta = "y", start = 0) +
    xlim(0.5, 2.5) + 
    scale_fill_manual(values = c(
      "Positive" = color_pos,
      "Negative" = "#E8E8E8"
    )) +
    annotate("text", x = 0.5, y = 0,
             label = paste0(mean_val, "%\n± ", sem_val),
             size = 4.5, fontface = "bold", color = "black") +
    labs(title = group_name, fill = "SA-beta-gal") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 13,
                                face = "bold", margin = margin(b = 5)),
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 9)
    )
}

p_pl  <- make_donut("PL",  "#A9A9A9")
p_sen <- make_donut("SEN", "#D73027")

p_pl + p_sen + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave('figures/Fig1d_RTcellline_SA-B-GAL_donutplot.pdf', width = 6, height=4, dpi=300)






