library(data.table)
library(cancereffectsizeR)
library(ggplot2)

# Load helper function to annotate noncoding variants.
helpers = new.env()
source('analysis/helpers.R', local = helpers)

cesa = load_cesa('output/cca_cesa.rds')

# Make top gene effects plot
rec_variants = cesa$variants[maf_prevalence > 1, variant_id]

tmp = variant_counts(cesa, rec_variants, by = 'study')[, -"variant_type"][, .SD, 
                                                                          .SDcols = patterns('variant_id|prevalence')]
melted = melt(tmp, id.vars = 'variant_id')[variable != 'total_prevalence']
multi_study_variants = melted[value > 0, uniqueN(variable), by = 'variant_id'][V1 > 1, variant_id]

effects = cesa$selection$all_effects[multi_study_variants, on = 'variant_id']

# Manuscript claim: 302 multi-study variants
stopifnot(effects[, .N == 302])

# Get descriptive labels for noncoding variants
to_anno = effects[is.na(gene), variant_id]
annotated = helpers$annotate_noncoding(cesa, to_anno, ces.refset.hg38$transcripts)
annotated[cesa$variants, display_name := paste0(chr, ':', 
                                                format(start, big.mark = ',', trim = TRUE),
                                                ' ', ref, '>', alt), on = 'variant_id']
annotated[gene == 'intergenic', variant_label := display_name]
annotated[gene != 'intergenic', variant_label := paste0(display_name, ' (', bare_anno, ')')]


effects[variant_type == 'aac', variant_label := sub('.* ', '', variant_name)]
effects[annotated, let(gene = i.gene, variant_label = i.variant_label), on = 'variant_id']
effects[cesa$variants, essential_splice := essential_splice, on = 'variant_id']

# Because of CES variant prioritization, everything that didn't get annotated already is splice-disrupting.
stopifnot(effects[variant_type == 'snv' & is.na(variant_label), all(essential_splice)])
effects[variant_type == 'snv' & is.na(variant_label), variant_label := paste0(variant_name, '\n(splice)')]


effects[cesa$variants, conseq := fcase(aa_ref == aa_alt & essential_splice == FALSE, 'Silent',
                                       essential_splice == TRUE, 'Splice-disrupting',
                                       aa_alt == 'STOP', 'Premature stop codon',
                                       aa_ref != aa_alt, 'Missense'), on = 'variant_id']

fp = effects[order(-selection_intensity)]
fp[is.na(conseq), conseq := 'Other']
fp[, conseq := factor(conseq, levels = c('Premature stop codon', 'Splice-disrupting', 'Missense', 'Other', 'Silent'))]
top_genes = fp[1:50, setdiff(gene, 'intergenic')]
fp2 = fp[, .SD[.I %in% 1:50 | gene %in% top_genes]]
fp2[gene == 'intergenic', gene := '(intergenic)']

gg = plot_effects(fp2, group_by = 'gene', topn = uniqueN(fp2$gene),
                  label_individual_variants = 'variant_label', color_by = 'conseq', 
                  legend_size_name = 'Variant prevalence', legend_color_name = 'Consequence',
                  viridis_option = 'C', x_title = 'Cancer effect', y_title = 'Gene',
                  legend_size_breaks = c(.01, .025, .05, .075), label_text_size = 2.5,
             legend.position = c(.88, .31), seed = 101011) + 
  guides(size = guide_legend(nrow = 2, byrow = T)) + 
  theme(legend.title = element_text(size = 9), legend.text = element_text(size = 8))
# ggsave('figures/top_effects_coding.pdf', gg, width = 1075 * 2.5, height = 700 * 2.5, 
#        dpi = 'retina', units = 'px', device = cairo_pdf)

saveRDS(gg, 'figures/variant_effects/top_effect_panel.rds')


# A claim for results section
fgfr2 = fp[gene == 'FGFR2', variant_id]
with_these_fgfr2 = samples_with(cesa, any_of = fgfr2)
cesa$samples[with_these_fgfr2, table(cca_type, fgfr2_fusion, exclude = NULL)]

#         fgfr2_fusion
# cca_type FALSE <NA>
#      IHC    12    3

