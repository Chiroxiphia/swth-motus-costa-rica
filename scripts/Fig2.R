# =========================================================
# FIGURE 2 — styled like original manuscript
# Same data, original colors, no mean diamond
# =========================================================

library(tidyverse)
library(readr)
library(ggplot2)

# ---------------------------------------------------------
# 1) Read Figure 2 data if needed
#    (skip this section if fig2_dat already exists in memory)
# ---------------------------------------------------------
# fig2_dat <- read_csv("Figure2_data_minimum_stopover_duration_spring2024.csv")

# ---------------------------------------------------------
# 2) Make sure site order is correct
# ---------------------------------------------------------
fig2_dat <- fig2_dat %>%
  mutate(
    tagging_site = factor(
      tagging_site,
      levels = c("Finca Cántaros", "Rancho Quemado")
    )
  )

# ---------------------------------------------------------
# 3) Labels with sample size for x axis
# ---------------------------------------------------------
site_counts <- fig2_dat %>%
  count(tagging_site) %>%
  mutate(label = paste0(tagging_site, "\n(n=", n, ")"))

x_labels <- setNames(site_counts$label, site_counts$tagging_site)

# ---------------------------------------------------------
# 4) Colors matching original-style figure
# ---------------------------------------------------------
site_colors <- c(
  "Finca Cántaros" = "#1b9e77",   # green/teal
  "Rancho Quemado" = "#2c7fb8"    # blue
)

# ---------------------------------------------------------
# 5) Build plot
#    - boxplots in manuscript-style colors
#    - jittered individual points
#    - NO mean diamond
# ---------------------------------------------------------
p_fig2 <- ggplot(fig2_dat, aes(x = tagging_site, y = min_stopover_days, color = tagging_site)) +
  geom_boxplot(
    width = 0.45,
    outlier.shape = NA,
    fill = NA,
    linewidth = 1
  ) +
  geom_jitter(
    width = 0.08,
    size = 2.8,
    alpha = 0.75
  ) +
  scale_color_manual(values = site_colors, guide = "none") +
  scale_x_discrete(labels = x_labels) +
  labs(
    x = NULL,
    y = "Minimum stopover duration (days)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    axis.line = element_line(linewidth = 1),
    axis.ticks = element_line(linewidth = 1)
  )

# ---------------------------------------------------------
# 6) Show plot in RStudio
# ---------------------------------------------------------
print(p_fig2)

# ---------------------------------------------------------
# 7) Save outputs
# ---------------------------------------------------------
ggsave(
  "Figure2_minimum_stopover_duration_spring2024_colored.png",
  plot = p_fig2,
  width = 6.5,
  height = 5,
  dpi = 600
)

ggsave(
  "Figure2_minimum_stopover_duration_spring2024_colored.pdf",
  plot = p_fig2,
  width = 6.5,
  height = 5
)

ggsave(
  "Figure2_minimum_stopover_duration_spring2024_colored.tiff",
  plot = p_fig2,
  width = 6.5,
  height = 5,
  dpi = 600,
  compression = "lzw"
)