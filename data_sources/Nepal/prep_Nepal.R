library(data.table)
library(readxl)

# Note: Fixed row 5862 in sheet 1 of original Excel file (had two values in one field)
nonsynonymous = as.data.table(read_excel('data_sources/Nepal/original/Nepaletal_Hepatology2017_suppl4,5.xlsx', 
                                         sheet = 1, skip = 1))
synonymous = as.data.table(read_excel('data_sources/Nepal/original/Nepaletal_Hepatology2017_suppl4,5.xlsx', 
                                         sheet = 2, skip = 1))
nepal = rbind(nonsynonymous, synonymous, fill = T)

# MAF contains data from other studies (Gao, Chan-ON, Zou). We got this data from original sources,
# so restrict to just new samples, which all have "Patient" instead of study name in their LibraryName.
nepal = nepal[LibraryName %like% 'Patient']
nepal = nepal[, .(Tumor_Sample_Barcode = paste0('nepal_', LibraryName),
                  Chromosome = Chr, Start_Position = Start, Reference_Allele = Ref, 
                  Tumor_Seq_Allele2 = Tumor_variant)]
nepal[, Chromosome := gsub('^chr', '', Chromosome)]



fwrite(nepal, 'data_sources/Nepal/nepal.maf.gz', sep = "\t")

# All Nepal samples are IHC fluke-negative.
# They screened for FGFR2 fusion (per page 21 of supplement PDF); negative for all WXS samples.
nepal_key = data.table(Unique_Patient_Identifier = unique(nepal$Tumor_Sample_Barcode), 
                       cca_type = 'IHC', fluke_status = 'negative', fgfr2_fusion = FALSE)
fwrite(nepal_key, 'study_sample_keys/nepal_sample_key.txt', sep = "\t")
