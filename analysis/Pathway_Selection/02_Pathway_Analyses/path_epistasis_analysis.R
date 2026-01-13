# Purpose: Run pairwise epistasis analysis on all pathways with effect estimates.

library(cancereffectsizeR)
library(data.table)
library(pbapply)
library(ggplot2)

cesa = load_cesa("output/cca_cesa.rds")
path_effects = fread('output/pathway/panCCA_redefined_pathway_effects.txt')
path_info = fread('output/pathway/redefined_pathway_info.txt')

pathway_defs = readRDS(file='output/pathway/redefined_pathway_defs.rds')

# Convert to a long-form data.table of variant-pathway associations
all_vars = lapply(names(pathway_defs), function(pid) {
  data.table(path_id = pid, variant_id = select_variants(cesa=cesa, variant_ids=pathway_defs[[pid]]))
})
all_vars_tb = rbindlist(all_vars)
# Remove variant_id prefix
setnames(all_vars_tb, sub("^variant_id\\.", "", colnames(all_vars_tb)))

# Define compound variants for each pathway
comp = define_compound_variants(cesa = cesa, variant_table = all_vars_tb, 
                                by = "path_id", merge_distance = Inf)

run_name = "path_epistasis"
# Run epistasis on all compound variants grouped by pathway
# 40 pathways -> ~15 mins
cesa_tmp = ces_epistasis(cesa = cesa, variants = comp, samples = cesa$samples, run_name = run_name)
# The epistasis results are stored in cesa$epistasis as a list
epistasis_dt = cesa_tmp$epistasis[[1]]
epistasis_dt[, variant_A := sub("\\.1$", "", variant_A)]
epistasis_dt[, variant_B := sub("\\.1$", "", variant_B)]
fwrite(epistasis_dt, file='./output/path_epi_dt.csv')

# Filter significant pathway pairs (e.g., p_epistasis < 0.05)
significant_epistasis = epistasis_dt[p_epistasis < 0.05]
significant_epistasis[, fdr := p.adjust(p_epistasis, method = 'fdr')]
significant_epistasis <- significant_epistasis[fdr < 0.05, ]

# Only get pathway pairs with higher prevalence
# fp = significant_epistasis[nAB/n_total > .02 | AB_epistatic_ratio < 1][order(AB_epistatic_ratio)]
# fp[, ab_label := format(round(AB_epistatic_ratio, 2))]

# cesR built-in plot:
# gg = plot_epistasis(fp, pairs_per_row = 5)[[1]] +
#   geom_text(aes(x = x, y = 5e4, label = ab_label), nudge_x = .3)
# fwrite(epistasis_dt,file='./output/epistasis_output/path_epi_panCCA.csv')

