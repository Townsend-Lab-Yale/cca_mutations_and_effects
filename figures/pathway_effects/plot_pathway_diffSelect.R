library(cancereffectsizeR)
library(data.table)
library(ggplot2)

# We'll mark cancer gene pathways
pw_info = fread('output/pathway/redefined_pathway_info.txt')
pw_info[is_cancer_gene_path == T, path_display_name := sub('\\((P\\d+)\\)$', '(\\1\u2605)', path_display_name)]


subtype_path_effects = readRDS('output/pathway/subtype_path_effects.rds')
nominal = rbindlist(lapply(subtype_path_effects,
                           function(x) {
                             x$testing[, fdr := p.adjust(p, method = 'fdr')]
                             return(x$testing[p < .05])
                           }), idcol = 'pair')
all_effects = fread('output/pathway/subtype_path_effects_with_diff_testing.txt')

fp = all_effects[nominal, on = c('pair', 'path_id')]
fp[pw_info, path_display_name := i.path_display_name, on = 'path_id']

# Leave out EHC when PHC or DCC are also significant
ehc_pw = fp[cca_type == 'EHC', path_id]
ehc_to_remove = fp[path_id %in% ehc_pw, any(c('ihc_phc', 'ihc_dcc') %in% pair), by = 'path_id'][V1 == T, path_id]
fp = fp[! (path_id %in% ehc_to_remove & pair == 'ihc_ehc')]

fp[, display_pair := fcase(pair %in%  c('ihc_ehc', 'ihc_dcc', 'ihc_phc'), 'iCCA vs. pCCA/dCCA',
                           pair == 'phc_dcc', 'pCCA vs. dCCA')]

#fp = unique(all_effects[best, on = c('pair', 'path_id')], by = c('cca_type', 'path_id'))
fp[, subtype_label := fcase(cca_type == 'IHC', 'Intrahepatic',
                            cca_type == 'EHC', 'Extrahepatic',
                            cca_type == 'PHC', 'Perihilar',
                            cca_type == 'DCC', 'Distal')]
fp[, subtype_color := fcase(cca_type == 'IHC', 'dodgerblue4',
                            cca_type == 'EHC', 'tan',
                            cca_type == 'PHC', 'darkolivegreen2', 
                            cca_type == 'DCC', 'gold1')]


fp[pw_info, let(pw_rank = pw_rank), on = 'path_id']
fp = fp[order(pw_rank)]
fp[, signif_label := fcase(fdr < .05, '*', default = '')]

for_signif = fp[signif_label != '', .(x = selection_intensity, path_display_name, pair, display_pair, signif_label)]
for_signif = for_signif[, .(x = max(fp$ci_high_95), signif_label = signif_label[1], display_pair = display_pair[1]), by = c('pair', 'path_display_name')]

fp = unique(fp, by = c('display_pair', 'path_id', 'cca_type'))

curr_signif = for_signif[display_pair %like% "iCCA"]
gg = plot_effects(fp, color_by = 'subtype_color', group_by = 'path_display_name', color_label = 'subtype_label', legend_color_name = 'Subtype',
                  x_title = 'Pathway cancer effects', y_title = 'Pathway', order_by_effect = FALSE,
                  label_individual_variants = FALSE) + 
  scale_size_continuous(name = 'Frequency of\npathway mutation', 
                        labels = scales::label_percent(accuracy = 1), 
                        range = c(1, 5),
                        guide = guide_legend(order = 2, ncol = 2, byrow = TRUE)) +
  facet_wrap(~display_pair, ncol = 2) + 
  theme(strip.text = element_text(size = 10),
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 11, family = 'Courier'),
        panel.spacing = unit(24, 'pt')) +
  guides(fill = guide_legend(order = 1, 
                             override.aes = list(size = 5))) 

curr_signif[gg@data, y := y_orig, on = 'path_display_name']


gg = gg + geom_text(data = curr_signif, aes(x = x, y = y, label = signif_label), 
            nudge_x = .2, nudge_y = .1, size = 14/.pt)

saveRDS(gg, 'figures/pathway_effects/path_diff_select_panel.rds')









