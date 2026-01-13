library(data.table)
library(readxl)
library(ces.refset.hg38)

# Read in Wardell Table S4
maf = as.data.table(read_excel('data_sources/Wardell/original/Wardell_mmc2-2.xlsx', sheet = 4, skip = 1))
maf = maf[, .(Unique_Patient_Identifier = paste0('wardell_', `Patient ID`), Chromosome = Chr, 
              Start_Position = Start_position, Reference_Allele, Tumor_Seq_Allele2)]

# Read in sample data. Sheet contains multiple tables; just taking the table for WGS/WXS samples.
sample_key = as.data.table(read_excel('data_sources/Wardell/original/Wardell_mmc2-2.xlsx', sheet = 1, 
                                      range = 'A2:M148'))
stopifnot(uniqueN(sample_key$ID) == sample_key[, .N])
sample_key = sample_key[, .(Unique_Patient_Identifier = paste0('wardell_', ID), `Sequencing method`, 
                            cca_type = Subtype, pM, viral = `Viral infection`, age = Age,
                            sex = Gender)]
sample_key[cca_type == 'ICC', cca_type := 'IHC']

# Leave out some gallbladder (GBC) and cystic duct (CDC) casees
sample_key = sample_key[cca_type %in% c('IHC', 'PHC', 'DCC')]
stopifnot(all(sample_key$pM == 0)) # all M0

# Wardell upplementary Figure 3 reports four patients with FGFR2 fusion. Only one of these is a WXS/WGS sample.
# Other fusion screening was seemingly performed but not reported.
sample_key[, fgfr2_fusion := Unique_Patient_Identifier == "wardell_RK279"]


# Table S4 includes calls from WXS, WGS, TGS. We're only using WXS/WGS from this source.
maf = maf[Unique_Patient_Identifier %in% sample_key$Unique_Patient_Identifier]

# Let's see if the reported somatic variants fall within the study's exome capture intervals.
wgs_maf = maf[Unique_Patient_Identifier %in% sample_key[`Sequencing method` == 'WGS', Unique_Patient_Identifier]]
wgs_check = preload_maf(maf = wgs_maf, refset = 'ces.refset.hg38', chain_file = 'reference/chains/hg19ToHg38.over.chain',
                        coverage_intervals_to_check = 'targeted_regions/Nextera_Rapid_Capture_hg38.bed.gz')
# > wgs_check[, mean(dist_to_coverage_intervals == 0)]
# [1] 0.9997139

wxs_maf = maf[! Unique_Patient_Identifier %in% wgs_maf$Unique_Patient_Identifier]
wxs_check = preload_maf(maf = wxs_maf, refset = 'ces.refset.hg38', chain_file = 'reference/chains/hg19ToHg38.over.chain',
                        coverage_intervals_to_check = 'targeted_regions/Nextera_Rapid_Capture_hg38.bed.gz')
# > wxs_check[, mean(dist_to_coverage_intervals == 0)]
# [1] 0.9998091

# It's clear that the reported calls for WGS samples have been trimmed to the intervals used in their exome sequencing.
# Therefore, we will treat all samples as exome-sequenced.
fwrite(maf, 'data_sources/Wardell/wardell.maf.gz', sep = "\t")



fwrite(sample_key[, -"Sequencing method"], 'study_sample_keys/wardell_sample_key.txt', sep = "\t")

