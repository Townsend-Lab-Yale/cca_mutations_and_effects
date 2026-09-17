library(data.table)

# All patient samples are in the sample key. And all samples are used in the CESAnalysis.
sk = fread('combined_sample_key.txt')
sk[Unique_Patient_Identifier %like% 'nakamura', study := 'nakamura']
sk[Unique_Patient_Identifier %like% 'jusakul', 
   study := paste0('jusakul.', coverage)]

# One sequencing type per study, with the above tweak
stopifnot(uniqueN(sk[, .(study, coverage)]) == uniqueN(sk$study))
counts = dcast(sk[, .N, by = c('study', 'cca_type')], 
               study ~ cca_type, fill = 0, value.var = 'N')
counts[, total := IHC + PHC + DCC + EHC]
setcolorder(counts, c('study', 'IHC', 'PHC', 'DCC', 'EHC', 'total'))

counts = counts[, s2 := tolower(study)][order(s2)][, s2 := NULL][]
counts = rbind(counts[study == 'yale'], counts[study != 'yale'])

fwrite(counts, 'output/study_included_sample_counts.txt', sep = "\t")

counts[, lapply(.SD, sum), .SDcols =  is.numeric]
# IHC   PHC   DCC   EHC total
# <int> <int> <int> <int> <int>
#   1:   930   186    74    15  1205

# More information for supplementary table
# 
kits_by_study = unique(sk[ , .(study, kit, coverage)])[, .(kits = .(paste(kit, collapse = ', '))), 
                                                       by = c('study', 'coverage')]
supp_info = copy(counts)
supp_info[kits_by_study, let(coverage = coverage, kits = kits), on = 'study']
substr(supp_info$study, 1, 1) = toupper(substr(supp_info$study, 1, 1))
supp_info[study == 'Jusakul.genome', let(study = 'Jusakul (WGS)',
                                         coverage = 'genome',
                                         kits = 'Not applicable')]
supp_info[study == 'Jusakul.targeted', let(study = 'Jusakul (TGS)',
                                           coverage = 'targeted',
                                           kits = 'Custom panel')]
setnames(supp_info, c('study', 'IHC', 'PHC', 'DCC', 'EHC', 'total'),
         c('source', 'iCCA_used', 'pCCA_used', 'dCCA_used', 'unspecified_eCCA_used',
           'num_samples_used'))

method_info = fread('data_sources/data_source_methods.txt')
supp_final = merge.data.table(supp_info, method_info, all = T, by = 'source')
setcolorder(supp_final, c('source', 'identifier', 'coverage', 'kits', 'num_samples_used'))
fwrite(supp_final, 'data_source_summary.txt', sep = "\t")
