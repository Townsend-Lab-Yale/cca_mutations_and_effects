library(data.table)

# Sorry, requires output files from figures/landscape/plot_landscapes.R.
msk_landscape = fread('output/landscape/msk_landscape_plotted_features.txt')
ihc_nontarget_landscape = fread('output/landscape/ihc_nontarget_landscape_plotted_features.txt')

get_summary_stats = function(dt) {
  melted = melt(dt[, .SD, .SDcols = patterns('(^pw)|(patient_id)')], id.vars = 'patient_id')
  melted[, patient_id := factor(patient_id, levels = unique(patient_id))]
  small_variants = melted[, .SD[value %like% '(SBS)|(Indel)', uniqueN(variable)], by = 'patient_id']
  sum1 = summary(small_variants$V1)
  all_variants = melted[, .SD[value != '', uniqueN(variable)], by = 'patient_id']
  sum2 = summary(all_variants$V1)
  
  tgs_only = summary(melted[, .SD[value != '' & variable %like% 'tgs', uniqueN(variable)], by = 'patient_id']$V1)
  exome_only = summary(melted[, .SD[value != '' & ! variable %like% 'tgs', uniqueN(variable)], by = 'patient_id']$V1)
  outputs = list(small_variants = sum1, all_variants = sum2, tgs_only = tgs_only, exome_only = exome_only)
  return(rbindlist(lapply(outputs, \(x) as.data.table(as.list(x))), idcol = 'variant_subset'))
}

msk_summary = get_summary_stats(msk_landscape)
ihc_nontarget_summary = get_summary_stats(ihc_nontarget_landscape)


output = rbindlist(list(MSK = msk_summary, IHC_nontarget = ihc_nontarget_summary), idcol = 'group')[order(variant_subset)]
fwrite(output, 'output/landscape/pathway_alteration_summary_stats.txt', sep = "\t")

# group variant_subset  Min. 1st Qu. Median     Mean 3rd Qu.  Max.
# <char>         <char> <num>   <num>  <num>    <num>   <num> <num>
# 1:           MSK   all_variants     0       1      2 2.767380       4    16
# 2: IHC_nontarget   all_variants     0       1      3 3.932755       5    34
# 3:           MSK     exome_only     0       0      0 0.000000       0     0
# 4: IHC_nontarget     exome_only     0       0      1 1.735358       2    19
# 5:           MSK small_variants     0       1      2 2.219251       3    16
# 6: IHC_nontarget small_variants     0       1      3 3.932755       5    34
# 7:           MSK       tgs_only     0       1      2 2.767380       4    16
# 8: IHC_nontarget       tgs_only     0       1      2 2.197397       3    17
