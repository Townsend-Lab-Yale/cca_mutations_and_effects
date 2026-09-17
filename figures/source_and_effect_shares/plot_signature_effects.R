library(cancereffectsizeR)
library(data.table)
library(patchwork)
library(cowplot)
library(ggtext)


source('figures/source_and_effect_shares/plot_signature_effects_helper.R') # plotting helper function

# Load in signature attributions
effect_list = readRDS(file = 'output/cca_signature_effects.rds') # from analyis/signature_effect_analysis.R
pan_effects = readRDS(file = 'output/cca_signature_effects_panCCA_only.rds')
effect_list[['panCCA']] = pan_effects



# Panels D and E were prepared in plot_variant_sources.R
other_stuff = readRDS('figures/source_and_effect_shares/for_signature_effects_figure.rds')
sig_info = other_stuff$sig_info
legend_info = other_stuff$legend_info
panel_d = other_stuff$panel_d
panel_e = other_stuff$panel_e


# Make source/effect share panels (Note: They exclude unknown-etiology and infrequently attributed signatures.)
scaled_width = 2 # panels must have same scale

gg_three_panel = plot_cca_signature_effects_three_panel(
  mutational_effects = effect_list,
  signature_groupings = sig_info,
  num_sig_groups = uniqueN(sig_info$description),
  other_color = "black",
  legend = plot_spacer(),
  panel_widths = c(scaled_width, -.02, scaled_width * 1.5, -.02, scaled_width * .5)
)

panel_widths = c(scaled_width, scaled_width * 1.5, scaled_width * .5)

# Need patchwork so that the scaling applies just to the plot areas
scaled_bars = Reduce(`+`, gg_three_panel) + plot_layout(nrow = 1, widths = panel_widths)


# Panel D: variant source probabilities for top variants
upper_panels = plot_grid(scaled_bars, plot_spacer() + theme_void(), rel_widths = c(.9, .05))

lower_panels = plot_grid(panel_d, panel_e, nrow = 2, label_size = 10,
         labels = c('D', 'E'), label_x = .006, label_y = c(.945, .965),
         rel_heights = c(.22, .78), scale = .95)

# Add top labels view cowplot since patchwork's tags wouldn't match up properly
full_figure_no_legend = plot_grid(upper_panels, lower_panels, ncol = 1, rel_heights = c(.47, .53)) +
  draw_label('A', .017, .933, size = 10, fontface = 'bold') +
  draw_label('B', .37, .933, size = 10, fontface = 'bold') + 
  draw_label('C', .78, .933, size = 10, fontface = 'bold')


# Homemade legend
ggs = list()
i = 1
for(grp in unique(legend_info$grp_category)) {
  fp = legend_info[grp_category == grp]
  # Thanks to https://stackoverflow.com/a/76835272 for legend sizing with ggtext
  legend_name = paste0("<span style='font-size: 9pt'>", grp, "</span>")
  if(i == 1) {
    legend_name = paste0("<span style='font-size: 12pt'><br><br>Signatures (COSMIC SBS)</span><br><br>", legend_name)
  }
  gg = ggplot(fp, aes(x = 0, fill = grp_color)) + geom_bar(color = 'black') +
    scale_fill_identity(drop = FALSE, guide = 'legend', labels = fp$full_label, 
                        breaks = fp$grp_color, name = legend_name) + 
    theme(legend.key.size = unit(.75, 'lines')) + theme_minimal()
  ggs[[i]] = gg + theme(legend.title = ggtext::element_markdown(),
                        legend.justification = 'top',
                        legend.box.just = 'top')
  i = i + 1
}

for_legend = Reduce(`+`, ggs) + plot_layout(guides = 'collect')
our_legend = plot_grid(get_legend(for_legend), plot_spacer() + theme_void(), ncol = 1, rel_heights = c(.82, .18))
full_figure = plot_grid(full_figure_no_legend, our_legend, 
                        nrow = 1, rel_widths = c(.74, .26), scale = c(1, .95))

ggsave('figures/prepped_figures/signature_effects.png', full_figure, width = 1325*2.5, height = 900*2.5, dpi = 'retina', 
       units = 'px', bg = 'white')



