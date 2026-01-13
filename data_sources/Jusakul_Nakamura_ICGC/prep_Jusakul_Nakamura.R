library(data.table)
library(readxl)
library(rtracklayer)
library(cancereffectsizeR)


jusakul_s1 = "data_sources/Jusakul_Nakamura_ICGC/original/Jusakul_S1.xlsx"
jusakul_s3 = "data_sources/Jusakul_Nakamura_ICGC/original/Jusakul_S3.xlsx"
jusakul_s4 = "data_sources/Jusakul_Nakamura_ICGC/original/Jusakul_S4.xlsx"
jusakul_s5 = "data_sources/Jusakul_Nakamura_ICGC/original/Jusakul_S5.xlsx"
nakamura_s2 = "data_sources/Jusakul_Nakamura_ICGC/original/Nakamura_S2.xlsx"

## Jusakul S1 describes WGS, WXS, and TGS samples
jusakul_key = as.data.table(read_excel(jusakul_s1, skip = 1, na = 'N/A'))[1:(.N - 2)] # leave out footer
setnames(jusakul_key, 'Sample ID', 'jusakul_id')

# We have 1 record per donor/sample
stopifnot(uniqueN(jusakul_key$jusakul_id) == jusakul_key[, .N])
jusakul_key$source = 'jusakul'
jusakul_key[endsWith(`TNM staging`, "M0"), pM := 0] 
jusakul_key[endsWith(`TNM staging`, "M1"), pM := 1]

jusakul_key[`Type (Fluke-Pos/Fluke-Neg)` == "Fluke-Neg", fluke_status := "negative"]
jusakul_key[`Type (Fluke-Pos/Fluke-Neg)` == "Fluke-Pos", fluke_status := "positive"]

jusakul_key[`Anatomical subtype` == "Intrahepatic", cca_type := "IHC"]
jusakul_key[`Anatomical subtype` == "Perihilar", cca_type := "PHC"]
jusakul_key[`Anatomical subtype` == "Distal", cca_type := "DCC"]
jusakul_key[`Anatomical subtype` == "Extrahepatic", cca_type := "EHC"] # a few unspecified EHC
jusakul_key[, Unique_Patient_Identifier := paste0('jusakul_', jusakul_id)]
setnames(jusakul_key, c('Age at surgery', 'Sex', 'Vital state (1=Dead)', 'Overall survival (days)'),
         c('age', 'sex', 'surv_status', 'surv_days'))

# The BTCA-SG MAF (all WGS) has some donors with multiple sequence samples. We will pick the samples that
# have identifiers corresponding to the Jusakul sample key. (In the MAF file, the sample identifiers have a "_T_" for tumor.)
jusakul_wgs_samples = jusakul_key[`Whole Genome sequencing` == "Yes", paste0('T_', jusakul_id)]
btca_sg_maf = fread("data_sources/Jusakul_Nakamura_ICGC/original/simple_somatic_mutation.open.BTCA-SG.tsv.gz")
stopifnot(uniqueN(btca_sg_maf[, .(icgc_sample_id, submitted_sample_id)]) == uniqueN(btca_sg_maf$icgc_sample_id))
jusa_wgs_maf = btca_sg_maf[jusakul_wgs_samples, on = 'submitted_sample_id', nomatch = NULL]

# We successfully found all 71 Jusakul WGS samples
stopifnot(uniqueN(jusa_wgs_maf$icgc_donor_id) == length(jusakul_wgs_samples),
          length(jusakul_wgs_samples) == 71)

# Merge the ICGC donor/sample IDs into the key file
jusa_wgs_maf[, jusakul_id := sub('^T_', '', submitted_sample_id)]
jusakul_key[jusa_wgs_maf, let(icgc_donor_id = icgc_donor_id, icgc_sample_id = icgc_sample_id), on = 'jusakul_id']

# ICGC's so-called MAF has duplicate records for all calls with multiple gene/transcript annotations, so
# we subset to unique calls.
jusa_wgs_maf = unique(jusa_wgs_maf[, .(Tumor_Sample_Barcode = paste0('jusakul_', jusakul_id), Chromosome = chromosome, Start_Position = chromosome_start, 
                                End_Position = chromosome_end, Reference_Allele = mutated_from_allele,
                                Tumor_Seq_Allele2 = mutated_to_allele)])
fwrite(jusa_wgs_maf, "data_sources/Jusakul_Nakamura_ICGC/btca_sg_wgs_jusakul.maf.gz", sep = "\t")


# Prepare Nakamura key and data. Nakamura and Jusakul keys will be combined at the end.
jusakul_key = jusakul_key[, .(Unique_Patient_Identifier, jusakul_id, icgc_donor_id, icgc_sample_id, fluke_status, cca_type, 
                              pM, HBV, HCV, Country, age, sex, surv_status, surv_days, TNM = `TNM staging`, source = 'jusakul')]
nakamura_key = as.data.table(read_excel(nakamura_s2, skip = 1))
stopifnot(uniqueN(nakamura_key$`Sample name`) == nakamura_key[, .N])
nakamura_key[, cca_type := fcase(`Tumor type 2` == 'Distal', 'DCC',
                                 `Tumor type 2` == 'IHCC', 'IHC',
                                 `Tumor type 2` == 'Perihilar', 'PHC', default = 'other')]

# A few samples are other cancer types (gallbladder cancer); exclude these.
nakamura_key = nakamura_key[cca_type != 'other']
nakamura_key[, `TNM staging` := paste0('T', nakamura_key$T, 'N', N, 'M', M)]
nakamura_key[`TNM staging` == "TNANNAMNA", `TNM staging` := NA] # convert pasted inclusions of "NA" to NAs
nakamura_key[endsWith(`TNM staging`, "M0"), pM := 0] 
nakamura_key[endsWith(`TNM staging`, "M1"), pM := 1]
setnames(nakamura_key, 'Sample name', 'nakamura_id')

# A few samples in the Nakamura key have a trailing T (like in the MAF)
nakamura_key[, nakamura_id := sub('T$', '', nakamura_id)]
nakamura_key[, Unique_Patient_Identifier := paste0('nakamura_', nakamura_id)]

# Get all available Nakamura CCA samples from the MAF.
# Just a few purportedly exome-sequenced samples from the study are not present; these were probably also
# analyzed with Jusakul WGS and therefore withheld from the WXS MAF.
btca_jp_wxs = fread('data_sources/Jusakul_Nakamura_ICGC/original/simple_somatic_mutation.open.BTCA-JP.tsv.gz')

# Just one sample per donor
stopifnot(uniqueN(btca_jp_wxs[, .(icgc_donor_id, icgc_sample_id)]) == uniqueN(btca_jp_wxs$icgc_donor_id))
btca_jp_wxs[, nakamura_id := sub('T$', '', submitted_sample_id)] # trim trailing T (for "tumor")
nakamura_wxs_maf = btca_jp_wxs[nakamura_key$nakamura_id, on = 'nakamura_id', nomatch = NULL]

# As with Jusakul, we subset to unique calls.
nakamura_wxs_maf = unique(nakamura_wxs_maf[, .(Tumor_Sample_Barcode = paste0("nakamura_", nakamura_id),
                                               Chromosome = chromosome, Start_Position = chromosome_start, 
                                               End_Position = chromosome_end, Reference_Allele = mutated_from_allele,
                                               Tumor_Seq_Allele2 = mutated_to_allele, icgc_donor_id, icgc_sample_id)])
nakamura_key[nakamura_wxs_maf, let(icgc_donor_id = icgc_donor_id, icgc_sample_id = icgc_sample_id), 
             on = c(Unique_Patient_Identifier = 'Tumor_Sample_Barcode')]

# Good news: No duplicate donors between Jusakul/Nakamura (at least, according to ICGC).
stopifnot(length(na.omit(intersect(nakamura_key$icgc_donor_id, jusakul_key$icgc_donor_id))) == 0)

# Check for unexpected overlap between Jusakul WGS and Nakamura WXS. The check_sample_overlap() function wants certain
# column names.
jusakul_for_check = copy(jusa_wgs_maf)
setnames(jusakul_for_check, c('Tumor_Sample_Barcode', 'Tumor_Seq_Allele2'), c('Unique_Patient_Identifier', 'Tumor_Allele'))
nakamura_for_check = copy(nakamura_wxs_maf)
setnames(nakamura_for_check, c('Tumor_Sample_Barcode', 'Tumor_Seq_Allele2'), c('Unique_Patient_Identifier', 'Tumor_Allele'))
poss_overlap = check_sample_overlap(maf_list = list(jusakul = jusakul_for_check, nakamura = nakamura_for_check))
poss_overlap = poss_overlap[variants_shared > 4] # see note below

# Apparently, 10 Jusakul WGS patients actually were part of the Nakamura WXS cohort, based on huge
# overlap (over 50% of Nakamura mutations found in corresponding Jusakul samples, for each of ten
# pairs; all pairs have over 50 shared variants). It is troubling that these samples have distinct
# ICGC IDs. Furthermore, Nakamura samples BD87 and BD89 appear to be from the same patient even
# though they have also have distinct ICGC IDs. Since BD87/BD89 do not have matching annotations (one
# listed as DCC, one as IHC), we will exclude both. The remaining pairs in poss_overlap are minor
# overlap that will be considered once the combined data set is assembled (see
# lift_filter_dedup_mafs.R).
stopifnot(poss_overlap[, .N] == 11)
nakamura_key = nakamura_key[! Unique_Patient_Identifier %in% c(poss_overlap$sample_A, poss_overlap$sample_B)]
nakamura_wxs_maf = nakamura_wxs_maf[nakamura_key$Unique_Patient_Identifier, on = "Tumor_Sample_Barcode", nomatch = NULL]
nakamura_key = nakamura_key[Unique_Patient_Identifier %in% nakamura_wxs_maf$Tumor_Sample_Barcode]

stopifnot(all(nakamura_key$`WES (Y=done, N-not done)` == 'Y')) # sanity check
nakamura_key = nakamura_key[, .(Unique_Patient_Identifier, nakamura_id, cca_type, pM, Virus, `HBV integration`, PSC, sex = Gender,
                                surv_days = `overall survival (days)`, TNM = `TNM staging`, source = 'nakamura')]
nakamura_key[surv_days == 'NA', surv_days := NA]
nakamura_key[, surv_days := as.numeric(surv_days)]

# Jusakul has some information on Nakamura samples that is not present in Nakamura!
# Therefore, we attempt to merge Nakamura/Jusakul sample data on shared columns (and we succeed).
nakamura_key_with_jusakul_info = merge.data.table(nakamura_key, jusakul_key[Country == 'Japan'], 
                                                  by = c('surv_days', 'cca_type', 'sex', 'TNM', 'pM'), 
                                                  allow.cartesian = TRUE)
nonunique_nakamura = nakamura_key_with_jusakul_info[, .N, by = 'nakamura_id'][N > 1, nakamura_id]

# We did well: 192/197 Nakamura samples have been uniquely paired with Jusakul identifiers.
paired_ids = nakamura_key_with_jusakul_info[! nonunique_nakamura, .(nakamura_id, jusakul_id), on = 'nakamura_id']
stopifnot(uniqueN(paired_ids) == 192)
nakamura_key[paired_ids, jusakul_id := jusakul_id, on = 'nakamura_id']

# Combine the sample keys
cols_to_merge = c(setdiff(names(jusakul_key), names(nakamura_key)), 'jusakul_id')
nakamura_key = merge.data.table(nakamura_key, jusakul_key[, .SD, .SDcols = cols_to_merge],
                                all.x = TRUE, all.y = FALSE, by = 'jusakul_id')
combined_key = rbind(nakamura_key, jusakul_key, fill = TRUE)

# Some of the patients are missing age, sex, or survival data. The UCSC Xena browser has some of this information
# for ICGC samples, so we'll merge it in.
icgc_survival = fread('data_sources/Jusakul_Nakamura_ICGC/ICGC_patient_info.tsv.gz')
setnames(icgc_survival, c('donor_age_at_diagnosis', 'donor_sex', 'donor_vital_status', 'donor_survival_time'),
         c('age_xena', 'sex_xena', 'surv_status_xena', 'surv_days_xena'))
icgc_survival = icgc_survival[icgc_donor_id %in% combined_key$icgc_donor_id]
icgc_survival[, sex_xena := fcase(sex_xena == 'female', 'F', sex_xena == 'male', 'M')]
icgc_survival[, surv_status_xena := fcase(surv_status_xena == 'alive', 0, surv_status_xena == 'deceased', 1)]
stopifnot(uniqueN(icgc_survival$icgc_donor_id) == icgc_survival[, .N])
combined_key = merge.data.table(combined_key, icgc_survival, all.x = TRUE, all.y = F, by = 'icgc_donor_id')

# sex_xena is consistent with sex. (Also, sex is never NA, so we don't need sex_xena.)
stopifnot(combined_key[sex_xena != sex, .N == 0], ! anyNA(combined_key$sex))
combined_key$sex_xena = NULL

## Age almost always matches. We'll set to NA where they conflict.
combined_key[age != age_xena, .(icgc_donor_id, age, age_xena)]
# icgc_donor_id   age age_xena
# 1:      DO217850    61       65
# 2:      DO218491    77       76
combined_key[age != age_xena, let(age = NA, age_xena = NA)]
combined_key[is.na(age), age := age_xena]
combined_key$age_xena = NULL

# When non-missing, vital status always matches. (Note that Nakamura didn't report vital status, so the matches are Jusakul.)
stopifnot(combined_key[! is.na(surv_status) & ! is.na(surv_status_xena), all(surv_status == surv_status_xena)])
combined_key[is.na(surv_status), surv_status := surv_status_xena]
combined_key$surv_status_xena = NULL

# Of patients with reports in both sources, 63 match and 5 don't
combined_key[, table(surv_days == surv_days_xena)] 

# Of the 5 patients with inconsistent survival reports, two are off by just a day; we'll keep the
# original reported days. For the other 3, we'll set all survival info to NA. (2 of these 3 also have NA age.)
combined_key[surv_days != surv_days_xena, .(age, surv_days, surv_days_xena)]
#       age surv_days surv_days_xena
# 1:    NA       149            176
# 2:    63       557            558
# 3:    45       212            578
# 4:    56      2099           2100
# 5:    NA       540            445

combined_key[abs(surv_days - surv_days_xena) > 1, 
             let(surv_days_xena = NA, surv_days = NA)]
combined_key[is.na(surv_days), surv_days := surv_days_xena]
combined_key$surv_days_xena = NULL

# Convert survival to months
days_in_month = 365.25/12
combined_key[, surv_month := surv_days / days_in_month]

# Already have HBV/HCV columns, and not using TNM in analysis. Country is encoded in jusakul_id and isn't needed.
combined_key[, c('Virus', 'HBV integration', 'TNM', 'Country', 'surv_days') := NULL] 

# Annotate fusions (Nakamura samples)
nakamura_fusions = as.data.table(readxl::read_excel("data_sources/Jusakul_Nakamura_ICGC/original/Nakamura_S2.xlsx",
                                                    sheet = "Supple. Table 7", skip = 1))
nakamura_fusions[, patient_id := paste0('nakamura_', sub('T$', '', Sample))]

combined_key[source == 'nakamura', let(fgfr2_fusion = FALSE, other_fusion = FALSE)]
nakamura_with_fgfr2 = nakamura_fusions[`Gene(5')` == 'FGFR2' | `Gene(3')` == "FGFR2", unique(patient_id)]
nakamura_with_other_fusion = nakamura_fusions[`Gene(5')` != 'FGFR2' & `Gene(3')` != "FGFR2", unique(patient_id)]
combined_key[nakamura_with_fgfr2, fgfr2_fusion := TRUE, on = 'Unique_Patient_Identifier']
combined_key[nakamura_with_other_fusion, other_fusion := TRUE, on = 'Unique_Patient_Identifier']

# Annotate fusions for Jusakul WGS (not reported for TGS samples)
jusakul_fusions = as.data.table(readxl::read_excel(jusakul_s4, skip = 1))
setnames(jusakul_fusions, c('Tumor ID', 'Gene at Left breakpoint', 'Gene at Right breakpoint'),
         c('Tumor', 'Gene1', 'Gene2'))
# Gene fusions are intra-chromosomal and chromosomal translocations that have gene annotations at each end
jusakul_fusions = jusakul_fusions[`SV type1` %in% c("ITX", "CTX") & ! is.na(Gene1) & ! is.na(Gene2)]
jusakul_fusions[, patient_id := paste0('jusakul_', Tumor)]
combined_key[source == 'jusakul' & Unique_Patient_Identifier %in% jusa_wgs_maf$Tumor_Sample_Barcode, 
             let(fgfr2_fusion = FALSE, other_fusion = FALSE)]
jusakul_with_fgfr2 = jusakul_fusions[Gene1 == 'FGFR2' | Gene2 == 'FGFR2', unique(patient_id)]
jusakul_with_other_fusion = jusakul_fusions[Gene1 != 'FGFR2' & Gene2 != 'FGFR2', unique(patient_id)]
combined_key[jusakul_with_fgfr2, fgfr2_fusion := TRUE, on = 'Unique_Patient_Identifier']
combined_key[jusakul_with_other_fusion, other_fusion := TRUE, on = 'Unique_Patient_Identifier']

# Nakamura samples underwent WXS using SureSelectV4 (+UTR) and SureSelect V5 (+UTR) panels. Here, we figure out which are which.
chain_file = "reference/chains/hg19ToHg38.over.chain" # UCSC chain files (licensed for non-commercial use only)

v4 = import.bed('targeted_regions/SureSelect_Human_All_Exon_V4_UTR_target_regions_hg38.bed.gz')
v5 = import.bed('targeted_regions/SureSelect_V5+UTR_hg38.bed.gz')
seqlevelsStyle(v5) = 'NCBI' # strip chr prefixes

# For clarity, remove calls in repetitive regions since they're more likely to be miscalls.
try_v4 = preload_maf(nakamura_wxs_maf, refset = 'ces.refset.hg38', coverage_intervals_to_check = v4, chain_file = chain_file)[is.na(problem)]
try_v4 = try_v4[repetitive_region == FALSE]
try_v5 = preload_maf(nakamura_wxs_maf, refset = 'ces.refset.hg38', coverage_intervals_to_check = v5, chain_file = chain_file)[is.na(problem)]
try_v5 = try_v5[repetitive_region == FALSE]

# Nakamura was mostly consistent about not making calls more than 100bp outside of target regions:
# For v4 and v5, check how many samples have 2 or more calls more than x bases outside covered regions
samples_with_calls = function(x, dt) {
  dt[dist_to_coverage_intervals > x, .N > 1,
     by = "Unique_Patient_Identifier"][V1 == T, uniqueN(Unique_Patient_Identifier)]
}
check = data.table(distance_outside = seq(from = 100, to = 1000, by = 50))
check[, with_v4 := sapply(distance_outside, samples_with_calls, try_v4)]
check[, with_v5 := sapply(distance_outside, samples_with_calls, try_v5)]

# > check
#     distance_outside with_v4 with_v5
# 1:              100     141      42
# 2:              150     135      14
# 3:              200     132       7
# 4:              250     129       6
# 5:              300     129       5
# 6:              350     127       4
# 7:              400     125       3
# 8:              450     124       3
# 9:              500     124       3
# 10:              550     123       3
# 11:              600     123       3
# 12:              650     122       3
# 13:              700     122       2
# 14:              750     122       2
# 15:              800     122       2
# 16:              850     122       1
# 17:              900     122       1
# 18:              950     121       1
# 19:             1000     120       1


# We will assume that the v5 samples are those with at least two calls >200bp out of coverage in v4
# regions. It's possible some low-TMB v5 samples will be misclassified as V4, but this is a pretty
# minimal source of error overall (especially since v4/v5 kits have high overlap).
v5_samples = try_v4[dist_to_coverage_intervals > 200, .N > 1, 
                    by = "Unique_Patient_Identifier"][V1 == T, unique(Unique_Patient_Identifier)]

v5_maf = nakamura_wxs_maf[Tumor_Sample_Barcode %in% v5_samples]
v4_maf = nakamura_wxs_maf[! Tumor_Sample_Barcode %in% v5_samples]

fwrite(v4_maf, 'data_sources/Jusakul_Nakamura_ICGC/nakamura_wxs_SureSelectV4-UTR.maf.gz')
fwrite(v5_maf, 'data_sources/Jusakul_Nakamura_ICGC/nakamura_wxs_SureSelectV5-UTR.maf.gz')

# Prepare Jusakul TGS data, ensuring exclusion of samples already present in Jusakul/Nakamura WGS/WXS.
## IMPORTANT: Silent mutations have been filtered out. The data is usable since the trinuc and gene
## mutation rates will be taken from WES/WGS for these samples, but the samples should be excluded for effect
## inferences on silent amino acid changes, if there are any of importance for the analysis.
tgs = as.data.table(read_excel(jusakul_s3, sheet = "TableS3A", skip = 1))
tgs = tgs[, .(Chromosome = Chr, Start_Position = `Start position`, Reference_Allele = Reference, Tumor_Allele = Variant,
              Tumor_Sample_Barcode = paste0("jusakul", "_", Tumor), jusakul_id = Tumor)]

# All but a few samples have Jusakul IDs.
jusakul_id_already_used = combined_key[unique(c(jusa_wgs_maf$Tumor_Sample_Barcode, nakamura_wxs_maf$Tumor_Sample_Barcode)), 
                                    jusakul_id, on = 'Unique_Patient_Identifier', nomatch = NULL]
tgs = tgs[! jusakul_id %in% jusakul_id_already_used]

# Now just need to make sure the 5 Nakamura samples without Jusakul IDs don't appear to match any Jusakul TGS samples.
samples_without_jusakul_id = combined_key[is.na(jusakul_id), Unique_Patient_Identifier]

n_for_check = nakamura_wxs_maf[samples_without_jusakul_id, on = 'Tumor_Sample_Barcode']
setnames(n_for_check, c('Tumor_Sample_Barcode', 'Tumor_Seq_Allele2'), c('Unique_Patient_Identifier', 'Tumor_Allele'))
tgs_for_check = copy(tgs)
setnames(tgs_for_check, 'Tumor_Sample_Barcode', 'Unique_Patient_Identifier')

# We'll do another sample overlap check with the full CCA combined data set, so for now just looking at
# overlap across the two groups.
poss_tgs_overlap = check_sample_overlap(maf_list = list(tgs = tgs_for_check, nk = n_for_check))
poss_tgs_overlap[source_A != source_B, .(sample_A, variants_A, sample_B, variants_B, variants_shared, greater_overlap)]
#          sample_A variants_A           sample_B variants_B variants_shared greater_overlap
# 1:   nakamura_BD5         74  jusakul_CCA_JP_15          3               3            1.00
# 2: nakamura_BD208        132 jusakul_CCA_JP_172          8               8            1.00
# 3: nakamura_BD210        206 jusakul_CCA_JP_174          8               6            0.75
# 4:   nakamura_BD5         74 jusakul_CCA_JP_180          2               1            0.50
# 5: nakamura_BD209         56 jusakul_CCA_JP_173          2               1            0.50

# Only 8 out of 177 Jusakul TGS samples have Japanese origin, yet all of our matches do.
# Out of an abundance of caution, we'll exclude all Japanese samples from TGS data.
# > unique(tgs_for_check[, .(jusakul_id)])[, table(jusakul_id %like% 'JP')]
# FALSE  TRUE 
# 169     8 
tgs = tgs[! jusakul_id %like% 'JP']
tgs$jusakul_id = NULL
fwrite(tgs, "data_sources/Jusakul_Nakamura_ICGC/jusakul_tgs.maf.gz", sep = "\t")

# Annotate that the TGS samples exclude synonymous mutations
combined_key[unique(tgs$Tumor_Sample_Barcode), excludes_synonymous := TRUE, on = "Unique_Patient_Identifier"]
all_used_samples = unique(c(tgs$Tumor_Sample_Barcode, nakamura_wxs_maf$Tumor_Sample_Barcode, jusa_wgs_maf$Tumor_Sample_Barcode))
stopifnot(! anyNA(all_used_samples))
combined_key = combined_key[Unique_Patient_Identifier %in% all_used_samples]

# Final check and save combined key
stopifnot(combined_key[, .N] == uniqueN(combined_key$Unique_Patient_Identifier),
          combined_key[, .N] == length(all_used_samples),
          ! anyNA(combined_key$cca_type))
fwrite(combined_key, 'study_sample_keys/jusakul_nakamura_sample_key.txt', sep = "\t")



# Build a targeted regions BED file for Jusakul TGS
jusakul_tgs_genes = as.data.table(read_excel(jusakul_s5, sheet = "TableS5C", skip = 1))

## Turns out 11 gene names in the Jusakul are not in Gencode; manually searched and found name conversions
jusakul_tgs_genes[, gencode_gene := Gene]
gencode_names  = c("GPAT4", "EMSY", "DOP1B", "ADGRG4", "ADGRV1", "H3C3",
                   "MRTFA", "KMT2D", "KMT2C", "CCN3", "PAK5")
names(gencode_names) = c("AGPAT6", "C11ORF30", "DOPEY2", "GPR112", "GPR98", "HIST1H3C", 
                         "MKL1", "MLL2", "MLL3", "NOV", "PAK7")
jusakul_tgs_genes[Gene %in% names(gencode_names), gencode_gene := gencode_names[Gene]]

# And 1 additional record not a gene, but rather a custom interval for the TERT promoter; set this aside for final BED output
tert_gr = GRanges(jusakul_tgs_genes[Gene == "TERT:PROMOTER", `Genomic region`])
jusakul_tgs_genes = jusakul_tgs_genes[Gene != "TERT:PROMOTER"]

stopifnot(all(jusakul_tgs_genes$gencode_gene %in% ces.refset.hg38$gene_names))

jusakul_tgs_genes[, chr := sub(':.*', '', `Genomic region`)]
jusakul_tgs_genes[, chr := sub('^chr', '', chr)]

# Via inspection of their target sizes, they appear to have captured CDS regions (as opposed to CDS + UTR)
gencode_exons = ces.refset.hg38$transcripts[type == 'CDS'][jusakul_tgs_genes, on = c('chr', gene_name = 'gencode_gene')]

# For each gene listed in Table S5, get all CDS regions and take intersection with their stated range
# (they just give a broad range from first base in their coverage through the last, rather than specific exons
# or targets)
chain_to_hg19 = import.chain('reference/chains/hg38ToHg19.over.chain')
names(chain_to_hg19) = sub('^chr', '', names(chain_to_hg19))
get_cds_grs = function(gene, jusakul_range) {
  gencode_gr = gencode_exons[gene_name == gene, reduce(unstrand(makeGRangesFromDataFrame(.SD)))]
  gencode_gr = reduce(unstrand(unlist(liftOver(gencode_gr, chain_to_hg19))))
  curr_gr = GRanges(jusakul_range)
  overlap_gr = reduce(intersect(gencode_gr, curr_gr))
  return(overlap_gr)
}
intersected_grs = mapply(get_cds_grs, jusakul_tgs_genes$gencode_gene, jusakul_tgs_genes$`Genomic region`)

# Turns out all their ranges had some overlap with Gencode CDS regions, and most were fairly close
jusakul_tgs_genes[, total_gencode_cds_size := sapply(intersected_grs, function(x) sum(width(x)))]

# Most of the time, our derived ranges are smaller than their stated "size", so our coverage ranges
# are conservative in terms of cancer effect inference.
jusakul_tgs_genes[, rel_size := total_gencode_cds_size / Size]
summary(jusakul_tgs_genes$rel_size)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.4402  0.8536  0.8849  0.8852  0.9242  1.2412 

# Combine all of our derived ranges, with the TERT promoter interval
final_grs = unlist(GRangesList(intersected_grs))
final_grs = c(final_grs, tert_gr)
seqlevelsStyle(final_grs) = "NCBI" # drop chr prefixes
final_grs = reduce(final_grs, drop.empty.ranges = T)

# See how well our BED intervals cover the TGS data.
cov_check = preload_maf(maf = tgs, refset = 'ces.refset.hg19', coverage_intervals_to_check = final_grs)
cov_check[, .(covered = mean(dist_to_coverage_intervals == 0),
              within_10 = mean(dist_to_coverage_intervals < 11),
              within_100 = mean(dist_to_coverage_intervals < 101))]
# covered within_10 within_100
# <num>     <num>      <num>
#   1: 0.9728162 0.9971004  0.9992751

## We'll supply 10bp padding to cover nearly all records
start(final_grs) = start(final_grs) - 10
end(final_grs) = end(final_grs) + 10
seqlevels(final_grs) = seqlevels(ces.refset.hg19$genome) # for sorting
final_grs = sort(final_grs)
rtracklayer::export.bed(final_grs, "targeted_regions/hg19/jusakul_tgs_hg19.bed.gz")

# Make the hg38 version, too.
lift_bed(bed = final_grs, chain = 'reference/chains/hg19ToHg38.over.chain', outfile = 'targeted_regions/jusakul_tgs_hg38.bed.gz')


# Gather focal CNA calls for some Jusakul samples, and for some Nakamura samples reported in Jusakul
# Jusakul reports 7 genes on one sheet (for WGS sample only) and one gene (ERBB2; most samples interrogated) on another
## Jusakul methods/supp say that this is relative copy number (i.e., changes relative to sample ploidy).
some_cna = as.data.table(readxl::read_excel(jusakul_s3, sheet = 'TableS3G', range = 'A2:H73'))
some_cna$Sample_ID = paste0("jusakul_", some_cna$Sample_ID) # convert identifiers to match project
setnames(some_cna, 'Sample_ID', 'Unique_Patient_Identifier')

some_cna = melt(some_cna, id.vars = 'Unique_Patient_Identifier', variable.name = 'gene', variable.factor = FALSE)
some_cna[, copy := fcase(value < 2, 'del', value == 2, 'neutral', value > 2, 'amp')]
some_cna = some_cna[, .(Unique_Patient_Identifier, gene, copy)]

# ERBB2 amplifications are reported separately. From Figure 1A in the paper, we can deduce that all
# other WGS samples had neutral copy state. We will interpret all TGS samples
# ("Validation" cohort) with ERBB2 not amplified as neutral state.
erbb2 = as.data.table(readxl::read_excel(jusakul_s3, sheet = 'TableS3E', skip = 2))
erbb2 = erbb2[, .(jusakul_id = sampleid, is_amp = `ERBB2 AMPLIFIED`)]

erbb2 = rbindlist(list(erbb2[is_amp == 'Yes', .(jusakul_id, gene = 'ERBB2', copy = 'amp')],
                       erbb2[is_amp == 'No', .(jusakul_id, gene = 'ERBB2', copy = 'neutral')],
                       erbb2[is_amp == 'NA', .(jusakul_id, gene = 'ERBB2', copy = NA)]))
erbb2[combined_key, Unique_Patient_Identifier := Unique_Patient_Identifier, on = 'jusakul_id']
erbb2 = erbb2[! is.na(Unique_Patient_Identifier)] # exclude samples not in key
erbb2$jusakul_id = NULL
all_cna = rbind(some_cna, erbb2)
all_cna = all_cna[! is.na(copy)]

stopifnot(all(all_cna$gene %in% ces.refset.hg38::ces.refset.hg38$gene_names))
stopifnot(all(all_cna$Unique_Patient_Identifier %in% combined_key$Unique_Patient_Identifier))
fwrite(all_cna, 'gene_copy_calls/Jusakul_Nakamura_gene_copy.txt', sep = "\t")

