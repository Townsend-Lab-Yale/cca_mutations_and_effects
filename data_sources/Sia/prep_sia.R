library(data.table)
maf = fread('data_sources/Sia/sia_CCA.maf.gz')

# Sia study was of IHC.
sample_key = data.table(Unique_Patient_Identifier = unique(maf$Tumor_Sample_Barcode), cca_type = 'IHC')


# Read in text version of Sia supplement's Table S1.
sia_fusions = fread('data_sources/Sia/Sia_TableS1_adapted.txt')
sia_fusions[, patient_id := paste0('sia_', Sample)]

with_fgr2 = sia_fusions[`5' gene` == 'FGFR2' | `3' gene` == 'FGFR2', unique(patient_id)]
with_other_fusion = sia_fusions[`5' gene` != 'FGFR2' & `3' gene` != 'FGFR2', unique(patient_id)]

sample_key[, fgfr2_fusion := Unique_Patient_Identifier %in% with_fgr2]
sample_key[, other_fusion := Unique_Patient_Identifier %in% with_other_fusion]
fwrite(sample_key, 'study_sample_keys/sia_sample_key.txt', sep = "\t")

