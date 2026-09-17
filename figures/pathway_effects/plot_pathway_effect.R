# Purpose: Plot pathway selection intensity for each CCA subtype and panCCA

library(cancereffectsizeR)
library(data.table)
library(ggplot2)

path_effects = fread('output/pathway/panCCA_redefined_pathway_effects.txt')[order(-selection_intensity)]
path_info = fread('output/pathway/redefined_pathway_info.txt')

# Add continuous color grouping column based on frequency
path_effects[, frequency := round(included_with_variant / (included_total + held_out), 3)]
path_effects[, frequency_pct := paste0(round(frequency * 100, 1), "%")]

# Original color: #DB382D
# Mark cancer gene pathways with stars (unicode \u2605)
path_info[is_cancer_gene_path == T, path_display_name := sub('\\((P\\d+)\\)$', '(\\1\u2605)', path_display_name)]
path_effects[path_info, path_display_name := path_display_name, on = 'path_id']
path_effects[, variant_name := path_display_name]

path_effects_modified <- copy(path_effects)
freq_range <- range(path_effects_modified$frequency, na.rm = TRUE)
freq_normalized <- (path_effects_modified$frequency - freq_range[1]) / diff(freq_range)

# Generate plasma colors for each frequency value
plasma_colors <- viridis::plasma(100)
color_indices <- pmax(1, pmin(100, round(freq_normalized * 99) + 1))
path_effects_modified[, point_fill_color := plasma_colors[color_indices]]

# Set all points to use the same color (this triggers the single legend behavior)
# But we'll override this in the scale_size guide
single_color <- plasma_colors[50]  # Use middle color as default

# Create the plot with manual color overrides in the size legend
gg <- plot_effects(path_effects_modified, 
                   group_by = 'path_display_name', 
                   x_title = 'Subpath cancer effect',
                   y_title = 'Subpath', 
                   legend.position = c(.72, .28), 
                   color_by = "frequency",
                   legend_size_name = 'Frequency of\n',
                   viridis_option = "plasma", 
                   topn = 40) +
  
  # Override the size scale to include color information
  scale_size_continuous(
    name = 'Frequency of\nsubpath mutation',
    labels = function(x) paste0(round(x * 100, 1), "%"),
    guide = guide_legend(
      title.position = 'top',
      ncol = 2,
      byrow = TRUE,
      override.aes = list(
        fill = viridis::plasma(6),
        alpha = 1
      )
    )
  ) +
  
  # Remove the fill legend since we're incorporating it into size
  guides(fill = "none") +
  
  theme(
    axis.title.x = element_text(size=10),
    axis.title.y = element_text(size=10),
    axis.text.x = element_text(size=8),
    axis.text.y = element_text(family = 'Courier', size=8),
    legend.title = element_text(size=10),
    legend.text = element_text(size=8),
  )
saveRDS(gg, 'figures/pathway_effects/path_effects_panel.rds')
