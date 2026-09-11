setwd('/data_hdd4/yuseong/cancer_aging/analysis/h.github/Fig5')

rm(list=ls())

library(dplyr)
library(openxlsx)
library(ggplot2)
library(patchwork)

### functions
build_curve <- function(conc, od) {
  std_avg <- aggregate(od, by = list(conc = conc), FUN = mean)
  std_avg <- std_avg[order(std_avg$x), ]  # sort by OD ascending for the inverse map
  splinefun(std_avg$x, std_avg$conc, method = "monoH.FC")
}

back_calc <- function(od, fit) {
  conc <- fit(od)
  pmax(conc, 0)
}

fmt_p <- function(p) if (p < 0.001) "P < 0.001" else sprintf("P = %.3f", p)

####


f <- "data/5l-o/MMP2_ELISA.xlsx"
std <- read.xlsx(f, sheet = "Standard_curve")
smp <- read.xlsx(f, sheet = "Samples")

# ---- Standard curve back-calculation ----
# A full 4-parameter logistic (4PL) fit was tried first, but the curve does not
# reach saturation within the tested range (OD is still rising steeply at the
# top 32 ng/mL standard), so the top asymptote is not identifiable and the fit
# is unstable (Top diverges to non-physical values, e.g. ~95, well above any
# observed OD). Since every unknown sample's OD falls well within the observed
# standard range, concentration is instead back-calculated by monotonic spline
# interpolation directly on the standard curve (mean OD of duplicates per
# point), which is standard, assumption-light practice when the curve top is
# not captured by the standard range.


experiments <- sort(unique(std$Experiment))
fits <- lapply(experiments, function(e) build_curve(std$Standard_ng_mL[std$Experiment == e],
                                                    std$OD[std$Experiment == e]))
names(fits) <- experiments

smp_avg <- smp %>%
  group_by(Experiment, Group, Cell_number, Media_volume_mL) %>%
  summarise(OD_mean = mean(OD), .groups = "drop") %>%
  rowwise() %>%
  mutate(Conc_ng_mL = back_calc(OD_mean, fits[[as.character(Experiment)]])) %>%
  ungroup()

blank_conc <- smp_avg %>% filter(Group == "Blank") %>% select(Experiment, Blank_conc = Conc_ng_mL)
smp_final <- smp_avg %>%
  filter(Group != "Blank") %>%
  left_join(blank_conc, by = "Experiment") %>%
  mutate(Conc_ng_mL_corr = pmax(Conc_ng_mL - Blank_conc, 0),
         Total_ng = Conc_ng_mL_corr * Media_volume_mL,
         pg_per_cell = ifelse(!is.na(Cell_number), Total_ng * 1000 / Cell_number, NA))

plot_df <- smp_final %>%
  mutate(
    pos = case_when(
      Group == "HUVEC_control"   ~ 1,
      Group == "HUVEC_senescent" ~ 2,
      Group == "FaDu_control"    ~ 3,
      Group == "FaDu_senescent"  ~ 4
    ),
    cell_type = ifelse(grepl("HUVEC", Group), "HUVEC", "FaDu")
  )

summary_df <- plot_df %>%
  group_by(pos, cell_type) %>%
  summarise(mean = mean(Conc_ng_mL_corr), se = sd(Conc_ng_mL_corr) / sqrt(n()), .groups = "drop")

points_df <- plot_df %>% transmute(pos, value = Conc_ng_mL_corr, cell_type)

p_huvec <- t.test(
  plot_df$Conc_ng_mL_corr[plot_df$Group == "HUVEC_control"],
  plot_df$Conc_ng_mL_corr[plot_df$Group == "HUVEC_senescent"],
  paired = TRUE
)$p.value

p_fadu <- t.test(
  plot_df$Conc_ng_mL_corr[plot_df$Group == "FaDu_control"],
  plot_df$Conc_ng_mL_corr[plot_df$Group == "FaDu_senescent"],
  paired = TRUE
)$p.value


ytop <- max(summary_df$mean + summary_df$se, points_df$value) * 1.45
y1 <- max(points_df$value[points_df$pos %in% c(1, 2)]) * 1.20
y2 <- max(points_df$value[points_df$pos %in% c(3, 4)]) * 1.20

main <- ggplot(summary_df, aes(x = pos, y = mean, fill = cell_type)) +
  geom_col(width = 0.6, color = "black") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.15) +
  geom_point(data = points_df, aes(x = pos, y = value), shape = 21,
             fill = "white", color = "black", size = 2, inherit.aes = FALSE) +
  annotate("segment", x = 1, xend = 2, y = y1, yend = y1, linewidth = 0.4) +
  annotate("text", x = 1.5, y = y1 * 1.10, label = fmt_p(p_huvec), size = 3.5) +
  annotate("segment", x = 3, xend = 4, y = y2, yend = y2, linewidth = 0.4) +
  annotate("text", x = 3.5, y = y2 * 1.10, label = fmt_p(p_fadu), size = 3.5) +
  scale_fill_manual(values = c("HUVEC" = "#9ab1bb", "FaDu" = "#e8b06c"), guide = "none") +
  scale_x_continuous(limits = c(0.4, 4.6), breaks = NULL, expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)), limits = c(0, ytop)) +
  labs(x = NULL, y = "Total MMP2 (ng/mL)") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text = element_text(color = "black"),
    plot.margin = margin(t = 10, r = 10, b = 2, l = 5)
  )

table_df <- data.frame(
  pos = c(1, 2, 3, 4, 1.5, 3.5),
  row = factor(c(rep("Senescent (IR)", 4), rep("Cell type", 2)),
               levels = c("Senescent (IR)", "Cell type")),
  label = c("-", "+", "-", "+", "HUVEC", "FaDu")
)

tbl <- ggplot(table_df, aes(x = pos, y = row, label = label)) +
  geom_text(size = 3.5) +
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

p <- main / tbl + plot_layout(heights = c(5, 1.3))
p

ggsave("figures/Fig5l_MMP2_ELISA_4cond.pdf", plot = p, width = 4, height = 6,
       units = "in", dpi = 300)
