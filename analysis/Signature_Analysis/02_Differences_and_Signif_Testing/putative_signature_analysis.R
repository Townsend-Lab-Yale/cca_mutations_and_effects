# Load in required libraries
library(data.table)
library(cancereffectsizeR)


# Load in bootstrapped attribution data (see signature_attribution_analysis.R):
# Contains the proportional attribution each signature contributes
# to each sample's total attribution, averaged across boots (370 samples)
signature_attributions = fread('output/final_unblended_signature_weights.txt')

sig_cols = names(signature_attributions)[names(signature_attributions) %like% 'SBS']

# Exclude signatures with no attributions (that is, artifact signatures and non-CCA signatures) from testing
all_zero_sigs = names(which(signature_attributions[, sapply(.SD, function(x) all(x == 0)), 
                                                   .SDcols = sig_cols]))
sig_cols = c(setdiff(sig_cols, all_zero_sigs), 'apobec', 'treatment')
signature_attributions[, apobec := SBS2 + SBS13]
signature_attributions[, treatment := SBS32 + SBS86 + SBS87 + SBS99]

signature_attributions = signature_attributions[, .SD, .SDcols = c(sig_cols, 'Unique_Patient_Identifier', 'cca_type')]

# Merge in signature detection frequency from bootstraps
raw_mp_out = fread('output/bootstrapped_mp_out.txt.gz')
raw_mp_out = raw_mp_out[Unique_Patient_Identifier %in% signature_attributions$Unique_Patient_Identifier]
raw_mp_out = raw_mp_out[, .SD, .SDcols = c('Unique_Patient_Identifier', 'boot',
                                           setdiff(sig_cols, c('treatment', 'apobec')))]
raw_mp_out[signature_attributions, cca_type := cca_type, on = 'Unique_Patient_Identifier']

melted = melt(raw_mp_out, id.vars = c('Unique_Patient_Identifier', 'boot', 'cca_type'), 
              variable.factor = FALSE, variable.name = 'signature')
putative = melted[, .(is_detected = mean(value > 0) >= .5), 
                  by = c('Unique_Patient_Identifier', 'signature')][is_detected == T, -"is_detected"]
putative_apobec = unique(putative[signature %in% c('SBS2', 'SBS13'),
                                  .(Unique_Patient_Identifier, signature = 'apobec')])
putative_treatment = unique(putative[signature %in% c('SBS32', 'SBS86', 'SBS87', 'SBS99'),
                                     .(Unique_Patient_Identifier, signature = 'treatment')])
putative = rbind(putative, putative_apobec, putative_treatment)

##### start of new code ####

apobec_samples = putative_apobec$Unique_Patient_Identifier # n = 179 samples with putative apobec signatures
treatment_samples = putative_treatment$Unique_Patient_Identifier # n = 178 samples with putative treatment signatures
haloalkane_samples = unique(putative[signature %in% c('SBS42')]$Unique_Patient_Identifier) # n = 134
aflatoxin_samples = unique(putative[signature %in% c('SBS24')]$Unique_Patient_Identifier) # n = 61
aristolochic_samples = unique(putative[signature %in% c('SBS22a', 'SBS22b')]$Unique_Patient_Identifier) # n = 39

# load in CESAnalysis object for more sample information
cesa = load_cesa('output/cca_cesa.rds')

treatment_table = table(cesa$samples[Unique_Patient_Identifier %in% treatment_samples]$study)
apobec_table = table(cesa$samples[Unique_Patient_Identifier %in% apobec_samples]$study)
haloalkane_table = table(cesa$samples[Unique_Patient_Identifier %in% haloalkane_samples]$study)
aflatoxin_table = table(cesa$samples[Unique_Patient_Identifier %in% aflatoxin_samples]$study)
aristolochic_table = table(cesa$samples[Unique_Patient_Identifier %in% aristolochic_samples]$study)

eligible_samples = cesa$samples[sig_analysis_eligible == TRUE, Unique_Patient_Identifier]
eligible_table = table(cesa$samples[Unique_Patient_Identifier %in% eligible_samples]$study)
all_names = names(eligible_table)

# Function to fill missing names with 0 and align
align_vector = function(vec, all_names) {
  out = setNames(rep(0, length(all_names)), all_names)
  out[names(vec)] = vec
  return(out)
}

# Align all vectors
treatment_vec = align_vector(treatment_table, all_names)
apobec_vec = align_vector(apobec_table, all_names)
haloalkane_vec = align_vector(haloalkane_table, all_names)
aflatoxin_vec = align_vector(aflatoxin_table, all_names)
aristolochic_vec = align_vector(aristolochic_table, all_names)
eligible_vec  = align_vector(eligible_table, all_names)  # should already be full

treatment_per_eligible = ifelse(eligible_vec == 0, NA, treatment_vec / eligible_vec)
apobec_per_eligible = ifelse(eligible_vec == 0, NA, apobec_vec / eligible_vec)
haloalkane_per_eligible = ifelse(eligible_vec == 0, NA, haloalkane_vec / eligible_vec)
aflatoxin_per_eligible = ifelse(eligible_vec == 0, NA, aflatoxin_vec / eligible_vec)
aristolochic_per_eligible = ifelse(eligible_vec == 0, NA, aristolochic_vec / eligible_vec)

# create output/putative_signature_group_detection.txt
sig_groups = list(
  apobec = apobec_vec,
  treatment = treatment_vec,
  haloalkane = haloalkane_vec,
  aflatoxin = aflatoxin_vec,
  aristolochic = aristolochic_vec
)

output_dt = data.table()
# loop through signature groups and studies
for (sig_name in names(sig_groups)) {
  sig_vec = sig_groups[[sig_name]]
  
  for (study_name in names(eligible_vec)) {
    output_dt = rbind(output_dt, data.table(
      sig_group = sig_name,
      study = study_name,
      detected_samples = sig_vec[study_name],
      total_samples = eligible_vec[study_name],
      detection_frequency = ifelse(eligible_vec[study_name] == 0, 
                                   NA, 
                                   sig_vec[study_name] / eligible_vec[study_name])
    ))
  }
  
  # add "all" rows that sum across studies
  output_dt = rbind(output_dt, data.table(
    sig_group = sig_name,
    study = "all",
    detected_samples = sum(sig_vec),
    total_samples = sum(eligible_vec),
    detection_frequency = sum(sig_vec) / sum(eligible_vec)
  ))
}

fwrite(output_dt, 'output/putative_signature_group_detection.txt', sep = '\t')

# run mutational_signature_effects
eligible_treatment_samples = intersect(eligible_samples, treatment_samples) # n = 178
eligible_apobec_samples = intersect(eligible_samples, apobec_samples) # n = 179
eligible_haloalkane_samples = intersect(eligible_samples, haloalkane_samples) # n = 134
eligible_aflatoxin_samples = intersect(eligible_samples, aflatoxin_samples) # n = 61
eligible_aristolochic_samples = intersect(eligible_samples, aristolochic_samples) # n = 39

sample_list = list(
  eligible_treatment_samples,
  eligible_apobec_samples,
  eligible_haloalkane_samples,
  eligible_aflatoxin_samples,
  eligible_aristolochic_samples
)

effect_list = list()
for (i in 1:length(sample_list)) {
  # calculate source shares and effect shares of mutational signatures.
  # iCCA
  iCCA_treatment_mut_effects = mutational_signature_effects(
    cesa,
    effects = cesa$selection$IHC,
    samples = intersect(sample_list[[i]], cesa$samples[cca_type == 'IHC', Unique_Patient_Identifier])
  )
  
  # pCCA
  pCCA_treatment_mut_effects = mutational_signature_effects(
    cesa,
    effects = cesa$selection$PHC,
    samples = intersect(sample_list[[i]], cesa$samples[cca_type == 'PHC', Unique_Patient_Identifier])
  )
  
  # dCCA
  dCCA_treatment_mut_effects = mutational_signature_effects(
    cesa,
    effects = cesa$selection$DCC,
    samples = intersect(sample_list[[i]], cesa$samples[cca_type == 'DCC', Unique_Patient_Identifier])
  )
  
  effects = list(iCCA = iCCA_treatment_mut_effects,
                 pCCA = pCCA_treatment_mut_effects,
                 dCCA = dCCA_treatment_mut_effects)
  
  effect_list = append(effect_list, effects)
}

treatment_effect_list = effect_list[1:3]
apobec_effect_list = effect_list[4:6]
haloalkane_effect_list = effect_list[7:9]
aflatoxin_effect_list = effect_list[10:12]
aristolochic_effect_list = effect_list[13:15]

# create output/putative_signature_group_shares.txt

signature_groups = list(
  apobec = c('SBS2', 'SBS13'),
  treatment = c('SBS32', 'SBS86', 'SBS87', 'SBS99'),
  haloalkane = c('SBS42'),
  aflatoxin = c('SBS24'),
  aristolochic = c('SBS22a', 'SBS22b')
)
cca_types = c('iCCA', 'pCCA', 'dCCA')
all_effects = list(
  apobec = apobec_effect_list,
  treatment = treatment_effect_list,
  haloalkane = haloalkane_effect_list,
  aflatoxin = aflatoxin_effect_list,
  aristolochic = aristolochic_effect_list
)

result_dt = data.table()
for (sig_group_name in names(signature_groups)) {
  for (cca_type in cca_types) {
    sigs = signature_groups[[sig_group_name]]
    effect_list = all_effects[[sig_group_name]]
    
    avg_source_share = sum(effect_list[[cca_type]]$mutational_sources$average_source_shares[sigs])
    avg_effect_share = sum(effect_list[[cca_type]]$effect_shares$average_effect_shares[sigs])
    result_dt = rbind(result_dt, 
                       data.table(sig_group = sig_group_name,
                                  cca_type = cca_type,
                                  average_source_share = avg_source_share,
                                  average_effect_share = avg_effect_share))
  }
}

# now run mutational_signature_effects() on all samples to get the signatures not represented

iCCA_mut_effects = mutational_signature_effects(
  cesa,
  effects = cesa$selection$IHC,
  samples = intersect(eligible_samples, cesa$samples[cca_type == 'IHC', Unique_Patient_Identifier])
)
pCCA_mut_effects = mutational_signature_effects(
  cesa,
  effects = cesa$selection$PHC,
  samples = intersect(eligible_samples, cesa$samples[cca_type == 'PHC', Unique_Patient_Identifier])
)
dCCA_mut_effects = mutational_signature_effects(
  cesa,
  effects = cesa$selection$DCC,
  samples = intersect(eligible_samples, cesa$samples[cca_type == 'DCC', Unique_Patient_Identifier])
)

# get remaining signatures
all_sigs = names(iCCA_mut_effects$mutational_sources$average_source_shares)
unrepresented_sigs = setdiff(all_sigs, unique(unlist(signature_groups)))

cca_mut_effects = list(
  iCCA = iCCA_mut_effects,
  pCCA = pCCA_mut_effects,
  dCCA = dCCA_mut_effects
)
for (sig in unrepresented_sigs) {
  for (cca_type in cca_types) {
    mut_effects = cca_mut_effects[[cca_type]]
    
    avg_source_share = mut_effects$mutational_sources$average_source_shares[sig]
    avg_effect_share = mut_effects$effect_shares$average_effect_shares[sig]
    
    result_dt = rbind(result_dt,
                       data.table(sig_group = sig,
                                  cca_type = cca_type,
                                  average_source_share = as.numeric(avg_source_share),
                                  average_effect_share = as.numeric(avg_effect_share)))
  }
}

fwrite(result_dt, 'output/putative_signature_group_shares.txt', sep = '\t')

