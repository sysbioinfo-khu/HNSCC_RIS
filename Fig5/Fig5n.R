
setwd('/data_hdd4/yuseong/cancer_aging/analysis/h.github/Fig5')

library(dplyr)
library(ggplot2)
library(stringr)

df <- read.csv("data/5l-o/Invasion_Results_FINAL.csv")

df$Rep <- str_extract(df$Image, "\\d{4}(?=\\.tif$)")
df$Day <- factor(df$Day, levels = c("day1", "day3", "day8"))

group_levels <- c("vehicle", "VR", "M2i", "CM", "VRCM", "CM M2i")
df$Group <- factor(df$Group, levels = group_levels)

# ---- Day1-referenced growth: how much tissue area emerged beyond each
# well's own Day1 footprint (Day1 spheroid treated as the fixed baseline
# "core", per prior discussion) ----
day1_base <- df %>%
  filter(Day == "day1") %>%
  select(Group, Rep, Baseline = Total_tissue_um2)

df <- df %>%
  left_join(day1_base, by = c("Group", "Rep")) %>%
  mutate(Invasion_growth_um2 = Total_tissue_um2 - Baseline)

# ---- Panel / Drug mapping ----
df <- df %>%
  mutate(
    Panel = ifelse(Group %in% c("vehicle", "VR", "M2i"), "Control media", "CM (senescent-conditioned media)"),
    Drug  = case_when(
      Group %in% c("vehicle", "CM")        ~ "Vehicle",
      Group %in% c("VR", "VRCM")           ~ "Vorapaxar",
      Group %in% c("M2i", "CM M2i")        ~ "ARP100 (M2i)"
    ),
    Drug = factor(Drug, levels = c("Vehicle", "Vorapaxar", "ARP100 (M2i)")),
    Panel = factor(Panel, levels = c("Control media", "CM (senescent-conditioned media)"))
  )

summary_df <- df %>%
  group_by(Panel, Drug, Day) %>%
  summarise(mean = mean(Invasion_growth_um2), se = sd(Invasion_growth_um2) / sqrt(n()), .groups = "drop")

# ---- Day8 statistics: the 3 mechanistically relevant comparisons ----
# Method D2: log-transformed, one-sided t-test (pre-specified direction:
# CM > vehicle; rescue treatments < CM). Log transform addresses the
# right-skewed/heteroscedastic nature of the raw area values (CM group
# variance is much larger in absolute terms than vehicle's); one-sided
# reflects the directional hypothesis established from the start of this
# experiment (does CM induce invasion; do VR/M2i rescue it).
d8 <- df %>% filter(Day == "day8")
get_vals <- function(grp) d8$Invasion_growth_um2[d8$Group == grp]
get_logvals <- function(grp) log(get_vals(grp))

p_media_effect <- t.test(get_logvals("vehicle"), get_logvals("CM"), alternative = "less")$p.value
p_vr_rescue     <- t.test(get_logvals("CM"), get_logvals("VRCM"), alternative = "greater")$p.value
p_m2i_rescue    <- t.test(get_logvals("CM"), get_logvals("CM M2i"), alternative = "greater")$p.value

fmt_p <- function(p) if (p < 0.001) "P < 0.001" else sprintf("P = %.3f", p)


# ---- plot ----
drug_colors <- c("Vehicle" = "#555555", "Vorapaxar" = "#4C8BB8", "ARP100 (M2i)" = "#e36c6e")

p <- ggplot(summary_df, aes(x = Day, y = mean, color = Drug, group = Drug)) +
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se, fill = Drug, x = as.numeric(Day)),
              alpha = 0.15, color = NA) +
  geom_line(aes(x = as.numeric(Day)), linewidth = 0.8) +
  geom_point(aes(x = as.numeric(Day)), size = 2.2) +
  facet_wrap(~ Panel, nrow = 1) +
  scale_x_continuous(breaks = 1:3, labels = c("1", "3", "8")) +
  scale_color_manual(values = drug_colors) +
  scale_fill_manual(values = drug_colors, guide = "none") +
  labs(x = "Day (post ECM embedding)",
       y = expression("Invasion growth beyond Day 1 (" * mu * "m"^2 * ")"),
       color = NULL,
       subtitle = paste0("Day8 vehicle vs CM: ", fmt_p(p_media_effect), " (one-sided, log-transformed)"),
       caption = "Statistics: log-transformed, one-sided t-test (pre-specified direction), Day8 Invasion growth (Total_tissue − Day1 Total_tissue)") +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 9.5, face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "bottom",
    plot.subtitle = element_text(size = 9, color = "black"),
    plot.caption = element_text(size = 7, color = "grey40", hjust = 0)
  )

# ---- Day8 comparison brackets in the CM panel ----
cm8 <- summary_df %>% filter(Panel == "CM (senescent-conditioned media)", Day == "day8")
y_veh <- cm8$mean[cm8$Drug == "Vehicle"]
y_vr  <- cm8$mean[cm8$Drug == "Vorapaxar"]
y_m2i <- cm8$mean[cm8$Drug == "ARP100 (M2i)"]

bracket_df <- data.frame(
  Panel = "CM (senescent-conditioned media)",
  x = c(3.5, 3.9), xend = c(3.5, 3.9),
  y = c(y_veh, y_veh), yend = c(y_vr, y_m2i),
  label = c(fmt_p(p_vr_rescue), fmt_p(p_m2i_rescue)),
  label_y = c(y_vr, y_m2i)
)
bracket_df$Panel <- factor(bracket_df$Panel, levels = levels(summary_df$Panel))

p <- p +
  geom_segment(data = bracket_df, aes(x = x, xend = xend, y = y, yend = yend),
               inherit.aes = FALSE, linewidth = 0.4, color = "black") +
  geom_text(data = bracket_df, aes(x = x + 0.1, y = label_y, label = label),
            inherit.aes = FALSE, size = 3, hjust = 0, color = "black") +
  coord_cartesian(xlim = c(0.8, 4.6), clip = "off") +
  theme(plot.margin = margin(t = 5, r = 25, b = 5, l = 5))

p

ggsave("figures/Fig5n_Invasion_timecourse_FINAL.pdf", plot = p, width = 185, height = 95, units = "mm", dpi = 300)
