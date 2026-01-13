library(cancereffectsizeR)
library(data.table)
library(ggplot2)

# Load in an annotation function (see helpers.R)
helpers = new.env()
source('analysis/helpers.R', local = helpers)


cesa = load_cesa('output/cca_cesa.rds')

subtype_variant_effects = readRDS('output/subtype_variant_effects.rds')
nominal = rbindlist(lapply(subtype_variant_effects,
                           function(x) {
                             x$testing[, fdr := p.adjust(p, method = 'fdr')]
                             return(x$testing[p < .05])
                           }), idcol = 'pair')

all_effects = fread('output/subtype_variant_effects_with_diff_testing.txt')
fp = all_effects[nominal, on = c('pair', 'variant_id')]

# Leave out EHC when PHC or DCC are also significant
ehc_variants = fp[cca_type == 'EHC', variant_id]
ehc_to_remove = fp[variant_id %in% ehc_variants, any(c('ihc_phc', 'ihc_dcc') %in% pair), by = 'variant_id'][V1 == T, variant_id]
fp = fp[! (variant_id %in% ehc_to_remove & pair == 'ihc_ehc')]

fp[, display_pair := fcase(pair %in%  c('ihc_ehc', 'ihc_dcc', 'ihc_phc'), 'iCCA vs. pCCA/dCCA/eCCA',
                           pair == 'phc_dcc', 'pCCA vs. dCCA')]

to_anno = fp[variant_type == 'snv', unique(variant_id)]
annotated = helpers$annotate_noncoding(cesa, to_anno, ces.refset.hg38$transcripts)
fp[variant_type == 'aac', variant_anno := variant_name]
fp[annotated, let(gene = i.gene, variant_anno = paste0(variant_name, '\n(', anno, ')')), on = 'variant_id']
fp[, best_gene_si := max(selection_intensity), by = 'gene']
fp = fp[order(-best_gene_si, -selection_intensity)]
fp[, subtype_label := fcase(cca_type == 'IHC', 'Intrahepatic',
                            cca_type == 'EHC', 'Extrahepatic',
                            cca_type == 'PHC', 'Perihilar',
                            cca_type == 'DCC', 'Distal')]
fp[, subtype_color := fcase(cca_type == 'IHC', 'dodgerblue4',
                            cca_type == 'EHC', 'tan',
                            cca_type == 'PHC', 'darkolivegreen2', 
                            cca_type == 'DCC', 'gold1')]

fp[, signif_label := fcase(fdr < .05, '*', default = '')]

# Remove LOF asterik in variant_anno
fp[, variant_anno := sub('\\*', "\u2731", variant_anno)] # heavy asterisk

# As stated in manuscript, excluding 15:61856541_A>C from plot, because it detected only in Nakamura, and 
# the sequencing context (A(10)C) suggests a likely calling error.
fp = fp[variant_id != '15:61856541_A>C']

# Create table that specifies where to put the significance markers on the plot.
for_signif = fp[signif_label != '', .(x = selection_intensity, variant_anno, pair, display_pair, signif_label)]
for_signif = for_signif[, .(x = max(fp$ci_high_95), signif_label = signif_label[1], display_pair = display_pair[1]), 
                        by = c('pair', 'variant_anno')]

fp = fp[order(-best_gene_si, gene, -(fdr < .05), -selection_intensity)]
fp = unique(fp, by = c('display_pair', 'variant_anno', 'cca_type'))


curr_signif = unique(for_signif[display_pair == "iCCA vs. pCCA/dCCA/eCCA"])

gg = plot_effects(fp[display_pair == "iCCA vs. pCCA/dCCA/eCCA"], group_by = 'variant_anno', 
                  color_by = 'subtype_color', color_label = 'subtype_label', legend_color_name = 'Subtype',
                  order_by_effect = FALSE, x_title = 'Cancer effects', y_title = 'Variant', 
                  label_individual_variants = FALSE) +
  scale_size_continuous(name = 'Variant prevalence', labels = scales::label_percent(accuracy = .1), range = c(1, 5)) +
  scale_x_log10(breaks = 10^(1:6), labels = c('  <10', format(10^(2:6), scientific = F, big.mark = ',', trim = T)),
               expand = expansion(add = c(.05, 0))) #+ coord_cartesian(expand = F, clip = 'off')

curr_signif[gg@data, y := y_orig, on = 'variant_anno']

gg = gg + geom_text(data = curr_signif, aes(x = x, y = y, label = signif_label), nudge_x = .2, nudge_y = .1, size=14/.pt) +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size=14),
        legend.title = element_text(size=14), 
        legend.text = element_text(size=14)) + 
  guides(fill = guide_legend(order = 1, 
                             override.aes = list(size = 5)),
         size = guide_legend(nrow = 2, byrow = T))

saveRDS(gg, 'figures/variant_effects/ihc-ehc_diff_select_panel.rds')

curr_signif = for_signif[display_pair == "pCCA vs. dCCA"] # empty, so no geom_text() this time
stopifnot(curr_signif[, .N] == 0)
# gg2 = plot_effects(fp[display_pair == "pCCA vs. dCCA"], group_by = 'variant_anno', 
#                   color_by = 'subtype_color', color_label = 'subtype_label', legend_color_name = 'Subtype',
#                   order_by_effect = FALSE, x_title = 'Cancer effects', y_title = 'Variant', 
#                   label_individual_variants = FALSE) +
#   scale_size_continuous(name = 'Variant frequency', labels = scales::label_percent(accuracy = 1), range = c(1, 5)) +
#   scale_x_log10(breaks = 10^(1:6), labels = c('  <10', format(10^(2:6), scientific = F, big.mark = ',', trim = T)),
#                 expand = expansion(0, 0)) + coord_cartesian(expand = F, clip = 'off') +
#   theme(text = element_text(size = 10),
#         axis.text.x = element_text(size = 10),
#         axis.text.y = element_text(size = 10),
#         axis.title.x = element_text(size = 10),
#         axis.title.y = element_text(size=10),
#         legend.title = element_text(size=10), 
#         legend.text = element_text(size=10)) + 
#   guides(fill = guide_legend(order = 1, 
#                              override.aes = list(size = 5)))
# 
# # Not currently saving these anywhere...
# # ggsave(gg2, file="figures/variant_diffSelect_phc-dcc.pdf", width = 4, height = 3, device = cairo_pdf)


