library(cancereffectsizeR)
library(data.table)
library(ggplot2)
library(cowplot)
library(extrafont)  

# May need to import Tahoma with extrafont::font_import() if systems fonts haven't previously been imported

cesa = load_cesa('output/cca_cesa.rds')

# Collect multi-study variants
rec_variants = cesa$variants[maf_prevalence > 1, variant_id]
study_counts = variant_counts(cesa, rec_variants, by = 'study')[, -"variant_type"][, .SD, 
                                                                          .SDcols = patterns('variant_id|prevalence')]
melted = melt(study_counts, id.vars = 'variant_id')[variable != 'total_prevalence']
multi_study_variants = melted[value > 0, uniqueN(variable), by = 'variant_id'][V1 > 1, variant_id]
effects = cesa$selection$all_effects[multi_study_variants, on = 'variant_id']

mut_effects = mutational_signature_effects(cesa, effects, samples = cesa$samples[sig_analysis_eligible == TRUE, Unique_Patient_Identifier])

prob_by_variant = copy(mut_effects$mutational_sources$average_by_variant)
prob_by_variant[effects, let(si = selection_intensity, 
                                                variant_name = variant_name), on = 'variant_id']
prob_by_variant = prob_by_variant[order(-si)][variant_id %in% multi_study_variants][1:20]



# Reproduce cannataro signature groupings and colors, with some tweaks
signature_groupings = list(
  "Deamination with age, clock-like" = "SBS1",
  "Unknown, clock-like" = "SBS5",
  "APOBEC" = c("SBS2", "SBS13"),
  "Defective homologous recombination" = "SBS3",
  "Tobacco" = c("SBS4", "SBS29"),
  "Chemotherapeutic agents" = c("SBS32", "SBS86", "SBS87", "SBS99"),
  "Alcohol-associated" = "SBS16",
  "Aristolochic acid exposure" = c("SBS22a", "SBS22b"),
  "Occupational haloalkane exposure" = "SBS42",
  "Aflatoxin exposure" = "SBS24",
  "Mismatch repair defects" = c("SBS15", "SBS20", "SBS21", "SBS44")
)

# Color mapping
cannataro_colors = c(
  "Deamination with age, clock-like" = "gray85",
  "Unknown, clock-like" = "gray60",
  "Infrequently attributed" = "gray35",
  "APOBEC" = "#7570b3",
  "Defective homologous recombination" = "#e7298a",
  "Tobacco" = "#a6761d",
  "Chemotherapeutic agents" = "#1b9e77",
  "Aflatoxin exposure" = "#579c9a",
  "Aristolochic acid exposure" = "#66a61e" ,
  "Occupational haloalkane exposure" = "#e6ab02",
  "Alcohol-associated" = "#d95f02",
  "Unknown etiology" = "black",
  "Mismatch repair defects" = "#8b324d"
)


non_exo_groups = c("Infrequently attributed", "Unknown etiology", "Deamination with age, clock-like",
                      "Unknown, clock-like", "APOBEC", "Mismatch repair defects", 
                      "Defective homologous recombination")
stopifnot(all(non_exo_groups %in% names(cannataro_colors)))

#stopifnot(all(names(cannataro_colors) %in% names(signature_groupings)))
samples_to_use = cesa$samples[sig_analysis_eligible == T, Unique_Patient_Identifier]
avg = cesa$mutational_signatures$biological_weights[samples_to_use, on = 'Unique_Patient_Identifier']
avg = avg[, lapply(.SD, mean), .SDcols = patterns('SBS')]

nonzero_sig = names(avg)[which(as.numeric(avg) != 0)]
nonzero_other = setdiff(nonzero_sig, unlist(signature_groupings))
all_unknown = cosmic_signature_info()[description == 'unknown etiology', name]

other_unknown = intersect(nonzero_other, all_unknown)
other_known = setdiff(nonzero_other, other_unknown)


avg[, variant_name := 'Average']
prob_by_variant[, c('variant_id', 'si') := NULL]
avg= avg[, .SD, .SD, .SDcols = names(prob_by_variant)]
prob_by_variant = rbind(avg, prob_by_variant)
prob_by_variant[, variant_name := factor(variant_name, levels = unique(variant_name))]

signature_groupings[['Infrequently attributed']] = other_known

for(i in 1:length(signature_groupings)) {
  curr_sbs = unname(signature_groupings[[i]])
  curr_grp_name = names(signature_groupings)[i]
  prob_by_variant[, (curr_grp_name) := rowSums(.SD), .SDcols = curr_sbs]
}
prob_by_variant = prob_by_variant[, .SD, 
                                  .SDcols = c('variant_name', names(signature_groupings))]


prob_by_variant[, `Unknown etiology` := 1 - rowSums(.SD), .SDcols = names(signature_groupings)]




melted = melt(prob_by_variant, id.vars = c('variant_name'), variable.factor = FALSE)
melted[, is_the_avg_row := variant_name == 'Average']

# Re-order colors to match multi-panel legend.
cannataro_colors = cannataro_colors[c("Chemotherapeutic agents", "Occupational haloalkane exposure", 
                                      "APOBEC", "Defective homologous recombination", "Aflatoxin exposure", 
                                      "Mismatch repair defects", "Aristolochic acid exposure", "Tobacco", 
                                      "Alcohol-associated", "Deamination with age, clock-like",
                                      "Unknown, clock-like", 
                                      "Infrequently attributed",
                                      "Unknown etiology")]
exo_groups = names(cannataro_colors)[! names(cannataro_colors) %in% non_exo_groups]
melted[, is_exo := variable %in% exo_groups]

plot_source_prob = function(dt, side = 'left', include_legend = FALSE) {
  color_levels = names(cannataro_colors)
  if(side != 'left') {
    color_levels = rev(names(cannataro_colors))
  }
  gg = ggplot(data = dt) + 
    geom_bar(aes(y = variant_name, weight = value, fill = variable, 
                 group = factor(variable, levels = color_levels)), 
             width = .8, color = 'black') + 
    scale_fill_manual(limits = names(cannataro_colors), values = cannataro_colors, labels = names(cannataro_colors), 
                      drop = FALSE) # guide = 'none
  if(side == 'left') {
    gg = gg + scale_x_reverse(expand = expansion(0, 0), limits = c(1.02, -0.02)) +
      labs(x = 'Variant-specific mutational source shares\n(Exogenous signatures)')
  } else if (side == 'right') {
    gg = gg + scale_x_continuous(expand = expansion(0, 0), limits = c(-0.02, 1.02)) +
      labs(x = 'Variant-specific mutational source shares\n(Endogenous and other signatures)')
  } else {
    stop('Bad side')
  }
  gg = gg + scale_y_discrete(limits = rev) +
    theme_minimal() + 
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.title.x = element_text(size = 8),
          plot.background = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())
  if(! include_legend) {
    gg = gg + theme(legend.position = 'none')
  }
  return(gg)
}


# Leaving out the known-other signatures because they're all small and sort of borderline as to whether
# they should be considered endogenous or exogenous
gg0a = plot_source_prob(melted[is_exo == TRUE & is_the_avg_row == TRUE], side = 'left') + 
  theme(axis.title.y = element_blank(), axis.title.x = element_text(size = 8)) + 
  labs(x = 'Average mutational source shares\n(Exogenous signatures)')
gg0b = plot_source_prob(melted[is_exo == FALSE & is_the_avg_row == TRUE], side = 'right') +
  theme(axis.title.y = element_blank(), axis.title.x = element_text(size = 8)) +
  labs(x = 'Average mutational source shares\n(Endogenous and other signatures)')
gg0_spacer = ggplot(data = data.table(x = 0, y = '')) +
  geom_text(aes(x, y, label = y), size = 2.75) +
  scale_x_continuous(expand = expansion(0, 0), breaks = 0, labels = '') +
  scale_y_discrete(limits = rev) +
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.title.x = element_text(size = 8),
        axis.text.y = element_blank(), axis.line = element_blank(),
        plot.subtitle = element_text(hjust = .7, size = 9, vjust = -2),
        panel.grid = element_blank(), axis.text.x = element_blank(),
        plot.background = element_blank(),
        plot.margin = unit(c(5.5, 0, 5.5, 0), units = 'pt')) + labs(x = '')

variants_by_exo_sum = melted[is_exo == T & variant_name != 'Average', sum(value), by = 'variant_name'][order(-V1), variant_name]
fp = melted[variants_by_exo_sum, on = 'variant_name'] # fp means for_plot :)

# Beautify SBS IDs
to_clean = cesa$variants[variant_type == 'snv'][unique(fp$variant_name), on = 'variant_name', nomatch = NULL]
to_clean[, for_display := paste0(chr, ':', format(start, big.mark = ','), ' ', ref, '>', alt)]

fp[to_clean, variant_name := for_display, on = 'variant_name']
fp[, variant_name := factor(variant_name, levels = unique(variant_name))]

gg1a = plot_source_prob(fp[is_exo == TRUE], side = 'left')
variant_labels = ggplot(data = fp[, .(x = 0, y = variant_name)]) + 
  geom_label(aes(x, y, label = y), size = 2.5, linewidth = 0, family = 'Tahoma') + 
  scale_x_continuous(expand = expansion(0, 0), breaks = 0, labels = '') +
  scale_y_discrete(limits = rev) +
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.title.x = element_text(size = 8),
        axis.text.x = element_text(size = 8), axis.text.y = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        plot.margin = unit(c(5.5, 0, 5.5, 0), units = 'pt')) + labs(x = '\n')
gg1b = plot_source_prob(fp[is_exo == FALSE], side = 'right')
  


important_widths = c(.48, .12, .48) # must be consistent for plot scales to match

panel_d = plot_grid(gg0a, gg0_spacer, gg0b, 
                     nrow = 1, rel_widths = important_widths)

panel_e = plot_grid(gg1a, variant_labels, gg1b, nrow = 1, rel_widths = important_widths)


# Build table of legend info
signature_groupings[['Unknown etiology']] = other_unknown
stopifnot(identical(sort(names(cannataro_colors)), sort(names(signature_groupings))))
signature_groupings = signature_groupings[names(cannataro_colors)]

sig_info = rbindlist(lapply(names(signature_groupings), function(desc) {
  sigs = signature_groupings[[desc]]
  data.table(
    name = sigs,
    short_name = as.character(sub("SBS", "", sigs)),
    description = desc,
    prioritize = FALSE,
    color = cannataro_colors[[desc]]
  )
}))
sig_info[, is_exo := description %in% exo_groups]



legend_info = unique(sig_info[, .(sig_group = description, grp_color = color, is_exo)])
legend_info[, grp_category := fcase(is_exo == TRUE, 
                                    'Environmental exposure signatures',
                                    default = 'Endogenous signatures')]
legend_info[sig_group == 'Unknown, clock-like', grp_category := 'Endogenous signatures']
legend_info[sig_group %in% c('Infrequently attributed', 'Unknown etiology'), grp_category := 'Other signatures']


sigs_by_group = sig_info[, for_label := paste0('(', paste0(short_name, collapse = ', '), ')'), by = 'description']



legend_info[sigs_by_group, for_label := for_label, on = c(sig_group = 'description')]

legend_info[, separator := fcase(nchar(paste0(sig_group, for_label)) > 45, '\n', default = ' ')]
legend_info[, full_label := paste0(sig_group, separator, for_label)]



# Save for inclusion in full figure
saveRDS(list(panel_d = panel_d, panel_e = panel_e, sig_info = sig_info, legend_info = legend_info), 
        'figures/source_and_effect_shares/for_signature_effects_figure.rds')




