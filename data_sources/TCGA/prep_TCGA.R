library(cancereffectsizeR)
library(ces.refset.hg38)


tcga_steele_ploidy = 'data_sources/TCGA/original/Steele S3 (41586_2022_4738_MOESM5_ESM).xlsx'

# Download samples from TCGA CCA primary tissue
# When last downloaded, acquired MAF data from Genomic Data Commons data release 40.0.
outname =  'data_sources/TCGA/TCGA-CHOL.maf.gz'
if(! file.exists(outname)) {
  get_TCGA_project_MAF(project = 'TCGA-CHOL', filename = outname, 
                       exclude_TCGA_nonprimary = TRUE)
}

maf = fread('data_sources/TCGA/TCGA-CHOL.maf.gz')

# Read clinical data, verifying that all MAF samples are present and no samples have conflicting information.
clinical = fread('data_sources/TCGA/clinical.project-tcga-chol.2024-08-01/clinical.tsv',
                 na.strings = "'--")
days_in_month = 365.25/12
clinical = clinical[, .(case_submitter_id, ajcc_pathologic_m, tissue_or_organ_of_origin, 
                        ajcc_pathologic_stage, sex = gender, surv_status = vital_status, 
                        age = age_at_diagnosis/365.25, surv_month = days_to_death / days_in_month)]

# Remove non-CCA samples
to_exclude = clinical[is.na(tissue_or_organ_of_origin) | tissue_or_organ_of_origin %in% c('Liver', 'Gallbladder'),
                      unique(case_submitter_id)]
maf = maf[! Unique_Patient_Identifier %in% to_exclude]
fwrite(maf, file = 'data_sources/TCGA/TCGA-CHOL.maf.gz', sep = "\t")

stopifnot(all(maf$Unique_Patient_Identifier %in% clinical$case_submitter_id))
stopifnot(uniqueN(clinical) == uniqueN(clinical$case_submitter_id))
clinical = unique(clinical)
clinical = clinical[case_submitter_id %in% maf$Unique_Patient_Identifier]

# summary(clinical$age)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   29.18   58.34   66.09   63.74   73.23   82.24       3 

clinical[, sex := fcase(sex == 'female', 'F', sex == 'male', 'M')]
# table(clinical$sex, exclude = NULL)
# F    M
#   25   21

clinical[, pM := ajcc_pathologic_m]
# When Stage is 1-3 and pM is NA or MX, set to M0
# > clinical[! pM %in% c('M0', 'M1'), table(ajcc_pathologic_stage)]
# ajcc_pathologic_stage
# Stage I  Stage II Stage IIB Stage III 
# 2         1         3         2 

clinical[pM == 'MX', pM := 'M0']
clinical[, pM := fcase(pM == 'M0', 0, pM == 'M1', 1)]
clinical[, c('ajcc_pathologic_m', 'ajcc_pathologic_stage') := NULL]
clinical[tissue_or_organ_of_origin == 'Intrahepatic bile duct', cca_type := 'IHC']
clinical[tissue_or_organ_of_origin == 'Extrahepatic bile duct', cca_type := 'EHC']
stopifnot(! anyNA(clinical$cca_type))
clinical$tissue_or_organ_of_origin = NULL
clinical[, surv_status := fcase(surv_status == 'Alive', 0, surv_status == 'Dead', 1)]
setnames(clinical, 'case_submitter_id', 'Unique_Patient_Identifier')

# TCGA portal doesn't have open-access fusion data, so we take it from the publication that introduced
# TCGA-CHOL, Farshidfar (10.1016/j.celrep.2017.02.033).
# Note: We checked GDC portal's protected STAR fusion data for these samples,
# and we confirmed that FGFR2 fusion status matches for all samples.
farshidfar_fusions = as.data.table(readxl::read_excel("data_sources/TCGA/original/Farshidfar_supp_tables.xlsx", 
                                                      sheet = "ST3. Fusions", skip = 1))

# Discard fusions that are marked "REJECT" (determined by Farshidfar to be unlikely to be functional)
farshidfar_fusions = farshidfar_fusions[verdict != 'REJECT']
farshidfar_fusions[, patient_id := substr(Sample, 1, 12)]
farshidfar_fusions = farshidfar_fusions[patient_id %in% clinical$Unique_Patient_Identifier]

tcga_with_fgfr = farshidfar_fusions[Gene1 == 'FGFR2' | Gene2 == "FGFR2", unique(patient_id)]
clinical[, fgfr2_fusion := Unique_Patient_Identifier %in% tcga_with_fgfr]

tcga_with_other_fusion = farshidfar_fusions[Gene1 != 'FGFR2' & Gene2 != "FGFR2", unique(patient_id)]
clinical[, other_fusion := Unique_Patient_Identifier %in% tcga_with_other_fusion]
fwrite(clinical, 'study_sample_keys/TCGA-CHOL_sample_key.txt', sep = "\t")


# From TCGA's GDC data portal. We'll compare integer copy numbers to ploidy.

file_to_id = fread('data_sources/TCGA/original/ASCAT3_gene_level_copy/gdc_sample_sheet.2024-08-19.tsv')
cna_files = list.files('data_sources/TCGA/original/ASCAT3_gene_level_copy', pattern = 'ascat3.*tsv.gz', full.names = TRUE)

file_to_id[, patient_id := substr(`Case ID`, 1, 12)]
filenames_to_match = sub('\\.gz', '', basename(cna_files))
names(cna_files) = file_to_id[filenames_to_match, patient_id, on = 'File Name', nomatch = NULL]

cna = rbindlist(lapply(cna_files, fread), idcol = 'Unique_Patient_Identifier')


# From Steele et al. (https://doi.org/10.1038/s41586-022-04738-6)
steele_ploidy_calls = as.data.table(readxl::read_excel(tcga_steele_ploidy, sheet = 'Genome doubling estimates'))
steele_ploidy_calls[, ploidy := fcase(GDx1 == 1 & GDx2 == 0, 4,
                                      GDx2 == 1, 2,
                                      default = 2)]
steele_ploidy_calls[, Unique_Patient_Identifier := substr(Sample, 1, 12)]
cna[steele_ploidy_calls, ploidy := ploidy, on = 'Unique_Patient_Identifier']

# The one sample not annotated in Steele is fully diploid
stopifnot(cna[is.na(ploidy), uniqueN(Unique_Patient_Identifier) == 1 & all(na.omit(copy_number == 2))])
cna[is.na(ploidy), ploidy := 2]

copies = cna[, .(Unique_Patient_Identifier, gene = gene_name, 
                 copy = fcase(copy_number > ploidy, 'amp', copy_number == ploidy, 'neutral', 
                              copy_number < ploidy, 'del'))][! is.na(copy)]

# We're not going to do any genes that have multiple entries per sample. Mostly snRNA genes, etc.
multiple_site_genes = copies[, .N, by = c('Unique_Patient_Identifier', 'gene')][N > 1, unique(gene)]
copies = copies[! gene %in% multiple_site_genes]

# Subset to coding genes in CES reference data, plus 3 additional genes used in the MSK data set (as
# seen prep_MSK.R).
copies = copies[gene %in% c(ces.refset.hg38$gene_names, "RYBP", "CDKN2B-AS1", "SOX2-OT")]

fwrite(copies, 'gene_copy_calls/TCGA-CHOL_gene_copy.txt', sep = "\t")
