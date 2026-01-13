library(data.table)

mut_effects = readRDS('output/cca_signature_effects.rds')
attributions = fread('output/final_unblended_signature_weights.txt')

all_signatures = setdiff(names(mut_effects$iCCA$effect_shares$by_sample), 'Unique_Patient_Identifier')

source_and_effect = melt(attributions, id.vars = c('Unique_Patient_Identifier', 'cca_type'), 
                         measure.vars = all_signatures, variable.name = 'signature', 
                         variable.factor = FALSE, value.name = 'source_share')

sample_effect_shares = rbindlist(lapply(mut_effects, function(x) x$effect_shares$by_sample), idcol = 'cca_type')
for_merge = melt(sample_effect_shares,  id.vars = c('Unique_Patient_Identifier', 'cca_type'), 
                 measure.vars = all_signatures, variable.name = 'signature', 
                 variable.factor = FALSE, value.name = 'effect_share')
for_merge[, cca_type := fcase(cca_type == 'iCCA', 'IHC',
                              cca_type == 'pCCA', 'PHC',
                              cca_type == 'dCCA', 'DCC', default = NA)]
stopifnot(! anyNA(for_merge$cca_type))

# Note that unspecified EHC samples will have NA effect_share. We're not using these samples since
# the effect shares are for iCCA, pCCA, dCCA.
source_and_effect[for_merge, effect_share := effect_share, on = c('Unique_Patient_Identifier', 'cca_type', 'signature')]

# Here's one test
treatment_signatures = c('SBS32', 'SBS86', 'SBS87', 'SBS99')
for_trt = source_and_effect[signature %in% treatment_signatures][,  
                                 .(cca_type = cca_type[1], 
                                   source_share = sum(source_share), 
                                   effect_share = sum(effect_share),
                                   signature = 'treatment'), by = 'Unique_Patient_Identifier']
apobec_signatures = c('SBS2', 'SBS13')
for_apobec = source_and_effect[signature %in% apobec_signatures, 
                               .(cca_type = cca_type[1], 
                                 source_share = sum(source_share), 
                                 effect_share = sum(effect_share),
                                 signature = 'apobec'), by = 'Unique_Patient_Identifier']

for_testing = rbindlist(list(source_and_effect, for_trt, for_apobec), use.names = TRUE)

for_testing = for_testing[cca_type != 'EHC'] # NA effect share on unspecified extrahepatic
subtype_testing = for_testing[, .(pval = wilcox.test(source_share, effect_share, paired = T, exact = F)$p.value),
                             by = c('cca_type', 'signature')]

pan_testing = for_testing[, .(cca_type = 'panCCA', pval = wilcox.test(source_share, effect_share, paired = T, exact = F)$p.value),
                              by = 'signature']
signif_testing = rbind(subtype_testing, pan_testing)
sorted_sbs = c('apobec', 'treatment', all_signatures)
signif_testing = signif_testing[sorted_sbs, on = 'signature']
signif_testing[, fdr := p.adjust(pval, 'fdr'), by = 'cca_type']

for_testing[signature %like% 'treat', mean(effect_share > source_share)]


# For manuscript
> signif_testing[signature %in% c('treatment', 'apobec', 'SBS16')]
# cca_type signature         pval          fdr
# <char>    <char>        <num>        <num>
# 1:      IHC    apobec 2.794187e-32 2.724332e-31
# 2:      PHC    apobec 1.764556e-11 1.676328e-10
# 3:      DCC    apobec 2.711135e-08 5.015600e-07
# 4:   panCCA    apobec 3.194910e-49 4.153384e-48
# 5:      IHC treatment 1.083062e-04 1.456531e-04
# 6:      PHC treatment 2.442762e-05 4.641248e-05
# 7:      DCC treatment 5.679583e-02 6.567018e-02
# 8:   panCCA treatment 1.063967e-08 1.595951e-08
# 9:      IHC     SBS16 4.185239e-01 4.579012e-01
# 10:      PHC     SBS16 1.454447e-05 2.908894e-05
# 11:      DCC     SBS16 9.303603e-06 2.458809e-05
# 12:   panCCA     SBS16 5.529538e-05 6.956516e-05

# Sanity check: Expecting effect share > source share
for_testing[signature  == 'treatment', mean(effect_share > source_share), by = 'cca_type']
# cca_type        V1
# <char>     <num>
# 1:      IHC 0.6104418
# 2:      PHC 0.6800000
# 3:      DCC 0.6341463


fwrite(signif_testing, 'output/source-vs-effect_signif.txt')


