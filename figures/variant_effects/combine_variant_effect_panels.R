library(ggplot2)
library(patchwork)

gg1 = readRDS('figures/variant_effects/top_effect_panel.rds') +
  theme(axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.text.x = element_text(hjust = .1),
        legend.position = 'right', legend.box.background = element_blank(),
        plot.tag.position = c(0, .87))

gg2 = readRDS('figures/variant_effects/ihc-ehc_diff_select_panel.rds') + 
  theme(legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        plot.tag.position = c(0, .84))

grid = gg1 + gg2 + plot_layout(guides = 'keep', ncol = 1, heights = c(.55, .45)) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 12, face = 'bold'),
        legend.position = 'right', legend.justification = 0)

ggsave(filename = 'figures/prepped_figures/variant_effects.png', plot = grid, width = 1000*2.5, height = 900*2.5, units = 'px')

