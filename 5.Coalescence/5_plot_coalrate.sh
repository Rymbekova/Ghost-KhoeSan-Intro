#!/usr/bin/env Rscript

library(relater)
library(ggplot2)
library(dplyr)
library(cowplot)

args = commandArgs(trailingOnly=TRUE)
dir = args[1]
prefix = args[2]

t_gen = 30

coal <- read.coal("relate_AA_all_chr_MaskHuman.coal") %>%
  group_by(group1, group2, epoch.start) %>%
  summarize(
    rate = mean(haploid.coalescence.rate, na.rm = TRUE),
    lower = quantile(haploid.coalescence.rate, 0.1, na.rm = TRUE),
    upper = quantile(haploid.coalescence.rate, 0.9, na.rm = TRUE),
    .groups = "drop"
  )

keep_groups <- c("admixing_group", "Pol_1", "Pol_2", "Pol_3", "Pol_4")
group_labels <- c(
  "admixing_group" = "Ghost",
  "Pol_1" = "San",
  "Pol_2" = "non-San",
  "Pol_3" = "Nea",
  "Pol_4" = "Den"
)

coal <- coal %>%
  filter(group1 %in% keep_groups, group2 %in% keep_groups) %>%
  mutate(
    group1 = group_labels[group1],
    group2 = group_labels[group2]
  )

# Left panel: within-group Ne (diagonal only), all lines overlapping
coal_diag <- coal %>% filter(group1 == group2)

p_diag <- ggplot(coal_diag) +
  geom_hline(yintercept = 1e4, linetype = "dashed", color = "grey") +
  geom_ribbon(aes(x = t_gen * epoch.start, ymin = 0.5 / upper, ymax = 0.5 / lower, fill = group1), alpha = 0.3, color = NA) +
  geom_step(aes(x = t_gen * epoch.start, y = 0.5 / rate, color = group1), linewidth = 1) +
  scale_x_continuous(trans = "log10") +
  annotation_logticks(sides = "bl") +
  scale_y_continuous(trans = "log10") +
  coord_cartesian(xlim = c(1e2, 1e6), ylim = c(1e2, 1e6)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  ) +
  labs(x = "years", y = "Ne", color = "Group", fill = "Group") +
  ggtitle("Effective population size within each group")

# Right panel: Ne to Ghost, all lines overlapping
coal_to_ghost <- coal %>%
  filter(group2 == "Ghost")

p_to_ghost <- ggplot(coal_to_ghost) +
  geom_hline(yintercept = 1e4, linetype = "dashed", color = "grey") +
  geom_ribbon(aes(x = t_gen * epoch.start, ymin = 0.5 / upper, ymax = 0.5 / lower, fill = group1), alpha = 0.3, color = NA) +
  geom_step(aes(x = t_gen * epoch.start, y = 0.5 / rate, color = group1), linewidth = 1) +
  scale_x_continuous(trans = "log10") +
  annotation_logticks(sides = "bl") +
  scale_y_continuous(trans = "log10") +
  coord_cartesian(xlim = c(1e2, 1e6), ylim = c(1e2, 1e6)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  ) +
  labs(x = "years", y = "Ne", color = "Group", fill = "Group") +
  ggtitle("Coalescence between each group and Ghost pop")

# Combine plots side-by-side
p_combined <- plot_grid(p_diag, p_to_ghost, ncol = 2, labels = c('A', 'B'), align = "h")

p_combined
ggsave(p_combined, file = paste0(dir, prefix, "relate_AA_all_chr_Ne_toGhost_overlap.pdf"), width = 10, height = 4.5)
