# Purpose: Plot epistatic interactions between pairs of genes across CCA types.
library(data.table)
library(stringr)
library(ggplot2)
library(ggrepel)
library(cancereffectsizeR)
library(cowplot)
library(grid)
library(gridExtra)

stopifnot(require(magick))

msk_exome_output_panCCA = fread(file='output/epistasis_output/gene_epi_MSKexome_panCCA.csv')
exome_output_panCCA = fread(file='output/epistasis_output/gene_epi_exome_panCCA.csv')

# variant A and B are always alphabetical (variant A < B).
# Excluding the exome inferences for genes with exome + MSK inferences.
all_epistasis = rbindlist(list(exome_tgs = msk_exome_output_panCCA,
                               exome_only = exome_output_panCCA[! msk_exome_output_panCCA, 
                                                                on = c('variant_A', 'variant_B')]), 
                          idcol = 'which_inference')
all_epistasis[, fdr := p.adjust(p_epistasis, method = 'fdr')]
fwrite(all_epistasis, 'output/epistasis_output/panCCA_epistatic_effects.txt', sep = "\t")


# Filter to FDR significance
signif_epistasis = all_epistasis[fdr < .05]
fp = signif_epistasis[nAB/n_total > .02 | AB_epistatic_ratio < 1][order(AB_epistatic_ratio)]

fp[, ab_label := format(round(AB_epistatic_ratio, 2))]

## Typical epistasis plot
# gg = plot_epistasis(fp, pairs_per_row = 5)[[1]] +
#   geom_text(aes(x = x, y = 5e4, label = ab_label), nudge_x = .3)
#   
# ggsave(file='figures/gene_epistasis_panCCA.png', gg,
#        width = 7, height = 5)

##### Epistasis plot alternative(s)
oncogene_lst <- c("KRAS", "NRAS", "IDH2")
tsuppressor_lst <- c("TP53", "CSMD3")
fp[, pair_id := .I]
fp[, driver_type_A := ifelse(variant_A %in% oncogene_lst, "oncogene",
                             ifelse(variant_A %in% tsuppressor_lst, "tumor suppressor",
                                    "non-driver gene"))]
fp[, driver_type_B := ifelse(variant_B %in% oncogene_lst, "oncogene",
                             ifelse(variant_B %in% tsuppressor_lst, "tumor suppressor",
                                    "non-driver gene"))]

# EPISTATIC TRAJECTORY PLOT
# Shows how each gene's effect changes from individual (A0/B0) to paired context (A_on_B/B_on_A)
create_epistatic_trajectory_plot <- function(data) {
  
  # Prepare data for trajectory plot using data.table
  dt_copy <- copy(data)
  
  # Calculate epistatic effect direction and magnitude for each pair
  dt_copy[, `:=`(
    delta_A = ces_A_on_B - ces_A0,
    delta_B = ces_B_on_A - ces_B0
  )]
  
  # Fill in missing data for lower bound of CI
  dt_copy[is.na(ci_low_95_ces_A_on_B), ci_low_95_ces_A_on_B := 0.001]
  dt_copy[is.na(ci_low_95_ces_B_on_A), ci_low_95_ces_B_on_A := 0.001]
  
  # Calculate FDR significance for coloring
  dt_copy[, fdr_significant := fdr < 0.05]
  
  # Create long format for trajectory plot
  # Gene A trajectory: ces_A0 -> ces_A_on_B
  dt_A_individual <- dt_copy[, .(
    pair_id, variant_A, variant_B, fdr_significant, ab_label,
    gene = variant_A,
    driver_type = driver_type_A,
    effect_size = ces_A0,
    ci_low = ci_low_95_ces_A0,
    ci_high = ci_high_95_ces_A0,
    context = "Individual",
    context_numeric = 1,
    p_change = p_A_change,
    delta = delta_A
  )]
  
  dt_A_paired <- dt_copy[, .(
    pair_id, variant_A, variant_B, fdr_significant, ab_label,
    gene = variant_A,
    driver_type = driver_type_A,
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
    pair_id, variant_A, variant_B, fdr_significant, ab_label,
    gene = variant_B,
    driver_type = driver_type_B,
    effect_size = ces_B0,
    ci_low = ci_low_95_ces_B0,
    ci_high = ci_high_95_ces_B0,
    context = "Individual",
    context_numeric = 1,
    p_change = p_B_change,
    delta = delta_B
  )]
  
  dt_B_paired <- dt_copy[, .(
    pair_id, variant_A, variant_B, fdr_significant, ab_label,
    gene = variant_B,
    driver_type = driver_type_B,
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
  
  # Handle zero or negative values for log scale by setting minimum threshold
  min_positive_value <- 0.01
  trajectory_data[effect_size <= 0, effect_size := min_positive_value]
  trajectory_data[ci_low <= 0, ci_low := min_positive_value]
  trajectory_data[ci_high <= 0, ci_high := min_positive_value]
  
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
  
  # Set the desired order for driver_type in the legend
  trajectory_data[, driver_type := factor(driver_type, 
                                          levels = c("oncogene", "tumor suppressor", "non-driver gene"))]
  
  # Determine which pairs go in each row (assuming ncol = 4)
  unique_pairs <- levels(trajectory_data$pair_id)
  n_pairs <- length(unique_pairs)
  ncol_facets <- 4
  
  # Pairs 1-4 go in row 1, pairs 5-7 go in row 2
  row1_pairs <- unique_pairs[1:min(4, n_pairs)]
  row2_pairs <- if(n_pairs > 4) unique_pairs[5:n_pairs] else character(0)
  print(row2_pairs)
  
  # Split data into two rows
  row1_data <- trajectory_data[pair_id %in% row1_pairs]
  row2_data <- if(length(row2_pairs) > 0) trajectory_data[pair_id %in% row2_pairs] else data.table()
  
  # Calculate unified y-axis limits across both rows for consistency
  all_values <- c(trajectory_data$effect_size, trajectory_data$ci_low, trajectory_data$ci_high)
  all_values <- all_values[all_values > 0]  # Remove any remaining non-positive values
  
  y_min <- min(all_values, na.rm = TRUE)
  y_max <- max(all_values, na.rm = TRUE)
  
  # Let ggplot2 handle log scale breaks and labels automatically
  
  # Function to create individual plot
  create_row_plot <- function(data, y_limits, n_col = 4, show_x_axis = TRUE, show_y_axis = TRUE, show_legend = TRUE) {
    
    # Only keep complete cases for plotting
    data_clean <- data[!is.na(effect_size) & 
                         !is.na(ci_low) & 
                         !is.na(ci_high) &
                         !is.na(context_numeric) &
                         effect_size > 0 &
                         ci_low > 0 &
                         ci_high > 0]
    
    if(nrow(data_clean) == 0) return(NULL)
    
    # Add x-axis jitter so genes A and B don't overlap
    data_clean[, x_jitter :=
                 context_numeric +
                 ifelse(gene == variant_A, -0.06,
                        ifelse(gene == variant_B, +0.06, 0))]
    
    # Create asterisk data for significant changes (only for paired context)
    asterisk_data <- data_clean[context == "Paired" & 
                                  !is.na(gene_change_significant) & 
                                  gene_change_significant == TRUE]
    
    # Prepare label data with jittering for overlapping cases
    label_data <- data_clean[context == "Paired"]
    
    # Calculate y-axis jitter based on log scale
    log_jitter_factor <- 2.5  # Multiplicative jitter for log scale
    
    # Apply y-axis jitter: when ab_label < 1, offset variant_A up and variant_B down
    label_data[, label_y_pos := effect_size]
    label_data[!is.na(ab_label) & ab_label < 1 & gene == variant_A, 
               label_y_pos := effect_size * log_jitter_factor]
    label_data[!is.na(ab_label) & ab_label < 1 & gene == variant_B, 
               label_y_pos := effect_size / log_jitter_factor]
    
    # Ensure label positions are within plot bounds
    label_data[label_y_pos > y_limits[2], label_y_pos := y_limits[2] / 1.1]
    label_data[label_y_pos < y_limits[1], label_y_pos := y_limits[1] * 1.1]
    
    # Create the plot
    p <- ggplot(data_clean, aes(x = x_jitter, y = effect_size)) +
      # Add confidence intervals
      geom_errorbar(aes(x = x_jitter, ymin = ci_low, ymax = ci_high,
                        color = driver_type),
                    width = 0.05, linewidth = 0.6, alpha = 0.5) +
      # Add lines connecting individual to paired effects
      geom_line(aes(x = x_jitter, group = interaction_group,
                    color = driver_type),
                linewidth = 1) +
      # Add filled points
      geom_point(aes(x = x_jitter, color = driver_type), 
                 size = 2.5, shape = 16) + 
      # Add asterisks for significant gene changes
      {if(nrow(asterisk_data) > 0) {
        geom_text(data = asterisk_data,
                  aes(x = x_jitter, y = effect_size),
                  label = "*", 
                  hjust = 1.75, vjust = 0.5, size = 14/.pt, fontface = "bold", 
                  color = "black", show.legend = FALSE)
      }} +
      # Add gene labels at the endpoints
      {if(nrow(label_data) > 0) {
        geom_text(data = label_data,
                  aes(x = x_jitter, y = label_y_pos, label = gene, color = driver_type), 
                  hjust = -0.2, size = 10/.pt, show.legend = FALSE)
      }} +
      # Faceting by gene pair
      facet_wrap(~pair_id, scales = "fixed", ncol = n_col) +
      scale_x_continuous(breaks = c(1, 2), 
                         labels = if(show_x_axis) c("Individual", "Paired") else c("", ""),
                         limits = c(0.8, 2.8)) +
      scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)),
                    limits = y_limits) +
      # Color palette
      scale_color_manual(values = c("oncogene" = "#E16A86", 
                                    "tumor suppressor" = "#009ADE", 
                                    "non-driver gene" = "darkgrey"),
                         labels = c("oncogene" = "oncogene",
                                    "tumor suppressor" = "tumor suppressor", 
                                    "non-driver gene" = "non-driver gene"),
                         name = "Gene classification") +
      labs(
        x = if(show_x_axis) "Somatic mutation context" else "",
        y = ""
      ) +
      theme_minimal() +
      theme(
        panel.grid.major.x = element_blank(),   # remove horizontal major gridlines
        panel.grid.minor.x = element_blank(),   # remove horizontal minor gridlines
        # keep vertical gridlines
        panel.grid.major.y = element_line(),
        panel.grid.minor.y = element_blank(),
        
        # keep axes
        axis.line = element_line(color = "grey", linewidth = 0.4),
        axis.ticks = element_line(color = "grey"),
        
        legend.position = if(show_legend) "right" else "none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 10),
        axis.title = element_blank(),
        strip.text = element_blank(),
        panel.spacing = unit(0.8, "lines"),
        plot.margin = if(!show_y_axis) margin(5.5, 5.5, 5.5, 2, "pt") else margin(5.5, 5.5, 5.5, 5.5, "pt")
      )
    
    return(p)
  }
  
  # Function to extract legend from a ggplot
  extract_legend <- function(plot) {
    tmp <- ggplotGrob(plot)
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }
  
  # Set unified y-axis limits for both rows
  unified_y_limits <- c(y_min * 0.8, y_max * 1.2)
  
  # Create row 1 plot (with legend to extract from)
  plot_row1_with_legend <- create_row_plot(
    data = row1_data,
    y_limits = unified_y_limits,
    n_col = 4,
    show_x_axis = FALSE,
    show_y_axis = FALSE,
    show_legend = TRUE
  )
  
  # Create row 1 plot without legend
  plot_row1 <- create_row_plot(
    data = row1_data,
    y_limits = unified_y_limits,
    n_col = 4,
    show_x_axis = FALSE,
    show_y_axis = FALSE,
    show_legend = FALSE
  )
  
  # Extract legend from the first plot
  legend <- extract_legend(plot_row1_with_legend)
  
  # Create row 2 plot (if data exists) - now using same scale as row 1
  plot_row2 <- if(nrow(row2_data) > 0) {
    create_row_plot(
      data = row2_data,
      y_limits = unified_y_limits,
      n_col = 3,  # Assume only 3 subplots in row 2
      show_x_axis = TRUE,
      show_y_axis = FALSE,
      show_legend = FALSE
    )
  } else {
    NULL
  }
  
  # Combine plots with custom layout
  if(!is.null(plot_row2)) {
    # Create layout matrix: 4 columns for row 1, 3 columns + 1 legend for row 2
    # Row 1: [plot1][plot1][plot1][plot1]
    # Row 2: [plot2][plot2][plot2][legend]
    
    # Convert plots to grobs
    grob_row1 <- ggplotGrob(plot_row1)
    grob_row2 <- ggplotGrob(plot_row2)
    
    # Create y-axis title as a separate text grob
    y_title <- textGrob("Scaled selection coefficient", 
                        rot = 90, 
                        gp = gpar(fontsize = 14),
                        x = 0.5, y = 0.5)
    
    # Create the main plot arrangement
    main_plots <- arrangeGrob(
      grob_row1, 
      arrangeGrob(grob_row2, legend, ncol = 2, widths = c(3, 1)),
      ncol = 1,
      heights = c(1, 1)
    )
    
    # Combine y-axis title with main plots
    combined_plot <- arrangeGrob(
      y_title, main_plots,
      ncol = 2,
      widths = c(0.05, 0.95),
      top = unit(0.5, "cm"),
      bottom = unit(0.5, "cm"),
      left = unit(0.5, "cm"),
      right = unit(0.5, "cm")
    )
    
    return(combined_plot)
  } else {
    return(plot_row1_with_legend)
  }
}

gene_epi_slope_plot <- create_epistatic_trajectory_plot(fp)
# ggsave(file='figures/gene_epistasis_panCCA_slopePlot.png', gene_epi_slope_plot,
#        width = 8, height = 5, dpi = 600)


# Read in network diagram (made in Illustrator)
network_panel = ggplot() + 
  cowplot::draw_image('figures/gene_epistasis/gene_epistasis_flow_diagram.png') +
  theme(panel.background = element_rect(fill = 'white'),
        plot.background = element_rect(fill = 'white'))

grid = cowplot::plot_grid(network_panel, gene_epi_slope_plot, ncol = 1,
                          labels = c('A', 'B'), label_y = c(.92, 1))
ggsave('figures/prepped_figures/gene_epistasis_plots.png', grid, height = 8, width = 8, bg = 'white')

