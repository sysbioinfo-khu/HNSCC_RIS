
setwd('/data_hdd4/yuseong/cancer_aging/analysis/h.github/Fig5')

rm(list=ls())

library(ggplot2)
library(dplyr)
library(ggsignif)

### functions
get_pval <- function(g1, g2) {
  v1 <- subset(df, Group == g1)$Value
  v2 <- subset(df, Group == g2)$Value
  if(sd(v1) == 0) {
    return(t.test(v2, mu = v1[1])$p.value)
  } else {
    return(t.test(v1, v2)$p.value)
  }
}

format_pval <- function(p) {
  if(p < 0.001) {
    return("P < 0.001")
  } else {
    return(sprintf("P = %.3f", p))
  }
}




df = read.csv('data/5l-o/permeability_result.csv')

df_summary <- df %>%
  group_by(Group) %>%
  summarise(
    Mean = mean(Value),
    SD = sd(Value),
    .groups = 'drop'
  )

treatment_labels <- c(
  "A" = "-\n-\n-",
  "B" = "-\n+\n-",
  "C" = "-\n-\n+",
  "D" = "+\n-\n-",
  "E" = "+\n+\n-",
  "F" = "+\n-\n+"
)

p_veh_sen <- get_pval("A", "D") # 1) Vehicle vs SEN CM + Vehicle
p_sen_vr  <- get_pval("D", "E") # 2) SEN CM + Vehicle vs SEN CM + VR
p_sen_m2i <- get_pval("D", "F") # 3) SEN CM + Vehicle vs SEN CM + M2i

anno_df <- data.frame(
  xmin = c(1, 4, 4),
  xmax = c(4, 5, 6),
  y_position = c(7.0, 5.4, 6.2), 
  annotations = c(format_pval(p_veh_sen), format_pval(p_sen_vr), format_pval(p_sen_m2i))
)


p <- ggplot() +
  geom_col(data = df_summary, aes(x = Group, y = Mean),
           fill = "#d3d3d3", color = "black", width = 0.6) +
  geom_errorbar(data = df_summary, aes(x = Group, ymin = Mean - SD, ymax = Mean + SD),
                width = 0.25, linewidth = 0.7) +
  geom_jitter(data = df, aes(x = Group, y = Value),
              width = 0.15, size = 2, alpha = 0.8, color = "black") +
  geom_signif(
    data = anno_df,
    aes(xmin = xmin, xmax = xmax, annotations = annotations, y_position = y_position),
    tip_length = 0.02, textsize = 4, vjust = -0.5, color = "black", manual = TRUE,
    inherit.aes = FALSE
  ) +
  scale_x_discrete(labels = treatment_labels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, 8.0)) +
  labs(
    y = "Relative Permeability\n(Fold change to Vehicle)",
    x = NULL
  ) +
  theme_classic() +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    axis.text.x = element_text(size = 12, lineheight = 1.5, color = "black", face = "bold"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 12)),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 110)
  )

p_final <- p + coord_cartesian(clip = "off") +
  annotate("text", x = 0.4, y = -0.3, 
           label = "SEN CM\nVorapaxar\nMMP2 inhibitor",
           hjust = 1, vjust = 1, size = 4.2, lineheight = 1.5, fontface = "bold", color = "black")

p_final

ggsave("figures/Fig5o_perm_near.pdf", plot = p_final, width = 6, height = 6,
       units = "in", dpi = 300)



