library(ggplot2)
library(ggrepel)
library(data.table)

# Pathway epistasis
epistasis_dt = fread(file='output/path_epi_dt.csv')

# Filter significant pathway pairs (e.g., p_epistasis < 0.05)
significant_epistasis = epistasis_dt[p_epistasis < 0.05]
significant_epistasis[, fdr := p.adjust(p_epistasis, method = 'fdr')]
significant_epistasis <- significant_epistasis[fdr < 0.05]

# Only get pathway pairs with higher prevalence

plot_epi_output = significant_epistasis[p_A_change < 0.05 | p_B_change < 0.05, ] # 3
setorder(significant_epistasis, AB_epistatic_ratio)
# plot_epistasis(plot_epi_output, variant_label_size = 6.5)

## Alternative plot
create_epistatic_trajectory_plot <- function(data) {
  
  # Prepare data for trajectory plot using data.table
  dt_copy <- copy(data)
  
  # Calculate epistatic effect direction and magnitude for each pair
  dt_copy[, `:=`(
    delta_A = ces_A_on_B - ces_A0,
    delta_B = ces_B_on_A - ces_B0
  )]
  
  # Fill in missing data for lower bound of CI
  dt_copy[is.na(ci_low_95_ces_A_on_B), ci_low_95_ces_A_on_B := .001]
  dt_copy[is.na(ci_low_95_ces_B_on_A), ci_low_95_ces_B_on_A := .001]
  
  # Calculate FDR significance for coloring
  dt_copy[, fdr_significant := fdr < 0.05]
  
  # Create long format for trajectory plot
  # Gene A trajectory: ces_A0 -> ces_A_on_B
  dt_A_individual <- dt_copy[, .(
    pair_id, variant_A, variant_B, fdr_significant,
    gene = variant_A,
    path_display_name = path_display_name_A,
    epi_display_name = epi_display_name_A,
    effect_size = ces_A0,
    ci_low = ci_low_95_ces_A0,
    ci_high = ci_high_95_ces_A0,
    context = "Individual",
    context_numeric = 1,
    p_change = p_A_change,
    delta = delta_A
  )]
  
  dt_A_paired <- dt_copy[, .(
    pair_id, variant_A, variant_B, fdr_significant,
    gene = variant_A,
    path_display_name = path_display_name_A,
    epi_display_name = epi_display_name_A,
    effect_size = ces_A_on_B,
    ci_low = ci_low_95_ces_A_on_B,
    ci_high = ci_high_95_ces_A_on_B,
    context = "Paired",
    context_numeric = 2,
    p_change = p_A_change,
    delta = delta_A
  )]
  
  # Gene B trajectory: ces_B0 -> ces_B_on_A
  dt_B_individual <- dt_copy[, .(
    pair_id, variant_A, variant_B, fdr_significant,
    gene = variant_B,
    path_display_name = path_display_name_B,
    epi_display_name = epi_display_name_B,
    effect_size = ces_B0,
    ci_low = ci_low_95_ces_B0,
    ci_high = ci_high_95_ces_B0,
    context = "Individual",
    context_numeric = 1,
    p_change = p_B_change,
    delta = delta_B
  )]
  
  dt_B_paired <- dt_copy[, .(
    pair_id, variant_A, variant_B, fdr_significant,
    gene = variant_B,
    path_display_name = path_display_name_B,
    epi_display_name = epi_display_name_B,
    effect_size = ces_B_on_A,
    ci_low = ci_low_95_ces_B_on_A,
    ci_high = ci_high_95_ces_B_on_A,
    context = "Paired",
    context_numeric = 2,
    p_change = p_B_change,
    delta = delta_B
  )]
  
  # Combine all trajectory data
  trajectory_data <- rbindlist(list(dt_A_individual, dt_A_paired, dt_B_individual, dt_B_paired), fill = TRUE)
  
  # Calculate max selection intensity by pair and create ordered factor for faceting
  trajectory_data[, max_selection := max(effect_size), by = pair_id]
  pair_order <- unique(trajectory_data[, .(pair_id, max_selection)])
  setorder(pair_order, -max_selection)
  # Convert pair_id to ordered factor based on max_selection (descending)
  trajectory_data[, pair_id := factor(pair_id, levels = pair_order$pair_id)]
  
  # Create interaction variable for grouping lines
  trajectory_data[, interaction_group := paste(pair_id, gene, sep = "_")]
  
  # Add significance for individual gene changes
  trajectory_data[, gene_change_significant := p_change < 0.05]
  
  # Only keep complete cases for plotting
  data_clean <- trajectory_data[!is.na(effect_size) & 
                                  !is.na(ci_low) & 
                                  !is.na(ci_high) &
                                  !is.na(context_numeric)]
  
  # Add x-axis jitter
  data_clean[, x_jitter :=
               context_numeric +
               ifelse(gene == variant_A, -0.06,
                      ifelse(gene == variant_B, +0.06, 0))]
  
  # Create a mapping from gene to path_display_name for the legend
  gene_to_path <- unique(data_clean[, .(gene, path_display_name)])
  gene_to_path <- gene_to_path[!is.na(gene) & !is.na(path_display_name)]
  
  # Set factor levels for gene based on unique values, and create labels mapping
  unique_genes <- unique(data_clean$gene)
  data_clean[, gene := factor(gene, levels = unique_genes)]
  
  # Create labels for legend (map from gene to path_display_name)
  legend_labels <- setNames(gene_to_path$path_display_name, gene_to_path$gene)
  
  # Create asterisk data for significant changes (only for paired context)
  asterisk_data <- data_clean[context == "Paired" & 
                                !is.na(gene_change_significant) & 
                                gene_change_significant == TRUE]
  
  # Prepare label data with jittering for overlapping cases
  label_data <- data_clean[context == "Paired"]
  label_data[, jitter_amount := 0.15]  # Adjusted for log scale
  
  # Apply jitter based on gene identity
  label_data[, label_y_pos := 10^((log10(effect_size) + log10(effect_size - delta))/2)]
  label_data[, label_x_pos := 1.5]
  # label_data[gene == variant_A, label_y_pos := effect_size * (1 + jitter_amount)]
  # label_data[gene == variant_B, label_y_pos := effect_size * (1 - jitter_amount)]
  
  # Generate colors for unique genes
  unique_genes_in_data <- unique(data_clean$gene)
  n_colors <- length(unique_genes_in_data)
  
  # Use a color palette that provides enough distinct colors
  if(n_colors <= 8) {
    colors <- RColorBrewer::brewer.pal(min(max(3, n_colors), 8), "Set2")[1:n_colors]
  } else {
    colors <- rainbow(n_colors)
  }
  names(colors) <- unique_genes_in_data
  
  # Create the plot
  gg <- ggplot(data_clean, aes(x = x_jitter, y = effect_size)) +
    # Add confidence intervals
    geom_errorbar(aes(x = x_jitter, ymin = ci_low, ymax = ci_high, color = gene),
                  width = 0.075, linewidth = 0.3, alpha = 0.75) +
    # Add lines connecting individual to paired effects
    geom_line(aes(x = x_jitter, group = interaction_group, color = gene), linewidth = 1) +
    # Add filled points
    geom_point(aes(x = x_jitter, color = gene), size = 3, shape = 16) +
    # Add asterisks for significant gene changes
    {if(nrow(asterisk_data) > 0) {
      geom_text(data = asterisk_data,
                aes(x = x_jitter, y = effect_size),
                label = "*", 
                hjust = 2, vjust = 0.6, size = 14/.pt, fontface = "bold", 
                color = "black", show.legend = FALSE)
    }} +
    # Add gene labels at the endpoints using epi_display_name
    {if(nrow(label_data) > 0) {
      geom_label_repel(data = label_data,
                       aes(x = label_x_pos, label = epi_display_name, color = 'black', y = label_y_pos), 
                       size = 2.5, show.legend = FALSE, alpha = .9)
    }} +
    # Faceting by gene pair
    facet_wrap(~pair_id, scales = "fixed", ncol = 3) +
    scale_x_continuous(breaks = c(1, 2), 
                       labels = c("Individual", "Paired"), 
                       limits = c(0.8, 2.8),
                       guide = guide_axis(n.dodge = 2)) +
    scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    # Color palette - use the legend_labels for the legend text
    scale_color_manual(values = colors,
                       labels = legend_labels,
                       name = "Subpath") +
    labs(
      x = "Somatic mutation context",
      y = "Scaled selection coefficient"
    ) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.position = "right",
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.text = element_text(family = 'Courier', size = 8),
      strip.text = element_blank(),
      panel.spacing = unit(0.8, "lines")
    )
  
  return(gg)
}

pw_info = fread('output/pathway/redefined_pathway_info.txt')
pw_info[, epi_display_name := paste0('P', pw_rank)]
pw_info[is_cancer_gene_path == TRUE, epi_display_name := sub('(^P\\d+)', '\\1\u2605', epi_display_name)]
plot_epi_output[pw_info, epi_display_name_A := i.epi_display_name, on=c('variant_A' = 'path_id')]
plot_epi_output[pw_info, epi_display_name_B := i.epi_display_name, on=c('variant_B' = 'path_id')]

pw_info[is_cancer_gene_path == TRUE, path_display_name := sub("(\\(P\\d+)\\)", "\\1\u2605)", path_display_name)]
plot_epi_output[pw_info, path_display_name_A := i.path_display_name, on=c('variant_A' = 'path_id')]
plot_epi_output[pw_info, path_display_name_B := i.path_display_name, on=c('variant_B' = 'path_id')]
# unique(c(plot_epi_output$path_display_name_A, plot_epi_output$path_display_name_B))

plot_epi_output[, pair_id := .I]
plot_epi_output[, ab_label := format(round(AB_epistatic_ratio, 2))]

gg = create_epistatic_trajectory_plot(plot_epi_output)
saveRDS(gg, 'figures/pathway_effects/path_epistasis_panel.rds')
# ggsave(file='figures/path_epistasis_panCCA.pdf', plot1,
#        width = 8, height = 3, device = cairo_pdf)
