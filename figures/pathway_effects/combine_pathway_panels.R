library(ggplot2)
library(patchwork)

path_diff = readRDS('figures/pathway_effects/path_diff_select_panel.rds')
path_effects = readRDS('figures/pathway_effects/path_effects_panel.rds')
path_epistasis = readRDS('figures/pathway_effects/path_epistasis_panel.rds')

path_diff = path_diff + 
  theme(title = element_blank(),
        axis.text.y = element_text(size = 10), 
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12), 
        legend.position = 'right',
        legend.title.position = 'top', legend.direction = 'vertical',
        legend.box.spacing = unit(0, 'pt'), legend.justification = 'left',
        legend.background = element_blank(), legend.box.background = element_blank())

path_effects = path_effects + 
  theme(title = element_blank(),
        axis.title.x = element_text(size = 12, margin = margin(9, 0, 0, 0, 'pt')),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))


path_epistasis = path_epistasis + 
  theme(title = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        axis.title.x = element_text(size = 12, margin = margin(9, 0, 0, 0, 'pt')),
        axis.title.y = element_text(size = 12))

upper_row = path_diff &
  plot_annotation(tag_levels = list('A')) & 
  theme(plot.tag.position = c(.01, .875))

lower_row = path_effects + path_epistasis +
  plot_layout(nrow = 1, guides = 'collect', widths = c(.41, .59)) &
  plot_annotation(tag_levels = list(c('B', 'C'))) &
  theme(plot.tag.position = c(.02, .93))

grid = upper_row + lower_row + 
  plot_layout(guides = 'keep', heights = c(.3, .7), ncol = 1) & 
  theme(plot.tag = element_text(size = 12, face = 'bold'),
        margins = margin())
ggsave('figures/prepped_figures/pathway_effects.png', grid, width = 8.8, height = 10)

