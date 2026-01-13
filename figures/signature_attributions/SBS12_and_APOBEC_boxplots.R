# Purpose: Visualize group comparisons of APOBEC and SBS12 signature attributions between CCA subtypes

# load in libraries
library(cancereffectsizeR)
library(data.table)
library(ggplot2)
library(cowplot)

# Load in bootstrapped attribution data (see main_analysis.R)
mp_out = fread('output/final_unblended_signature_weights.txt')

set.seed(8142023) # for reproducible jitter

plot_data = mp_out[cca_type %in% c("IHC", "PHC", "DCC")] # n = 365
plot_data$APOBEC = (plot_data$SBS2 + plot_data$SBS13)

# remap labels for proper CCA nomenclature
label_map = c("IHC" = "iCCA", "PHC" = "pCCA", "DCC" = "dCCA")
plot_data[, cca_display := factor(label_map[cca_type], levels = c("iCCA", "pCCA", "dCCA"))]


# Create the APOBEC plot
gg1 = ggplot(plot_data, aes(x = cca_display, y = APOBEC, fill = cca_display)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.3, alpha = 0.6, size = 1.5) +
  scale_fill_manual(values = c("iCCA" = "dodgerblue4", "pCCA" = "darkolivegreen2", "dCCA" = "gold1")) +
  theme_classic() +
  labs(x = "CCA subtype", y = "Proportion of tumor substitutions\nattributed to SBS2 + SBS13") +
  theme(legend.position = "none")

#ggsave(file = 'figures/APOBEC_signatures.png', plot = gg, width = 2700, height = 1300, units = 'px')

# SBS12 plot with same y-axis scaling
gg2 = ggplot(plot_data, aes(x = cca_display, y = SBS12, fill = cca_display)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.3, alpha = 0.6, size = 1.5) +
  scale_fill_manual(values = c("iCCA" = "dodgerblue4", "pCCA" = "darkolivegreen2", "dCCA" = "gold1")) +
  theme_classic() +
  labs(x = "CCA subtype", y = "Proportion of tumor substitutions\nattributed to SBS12") +
  theme(legend.position = "none")

grid = cowplot::plot_grid(gg1, gg2, ncol = 1, labels = c('A', 'B'), scale = .95)


ggsave(file = 'figures/prepped_figures/combined_APOBEC_SBS12_signatures.png', 
       plot = grid, width = 7, height = 8, units = 'in', bg = 'white')
