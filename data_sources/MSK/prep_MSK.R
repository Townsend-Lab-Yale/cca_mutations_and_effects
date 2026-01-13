library(data.table)
library(ces.refset.hg38)
stopifnot(packageVersion('ces.refset.hg38') >= as.package_version('1.3.0'))

# Note: MSK data set last downloaded 07-22-24 from https://cbioportal-datahub.s3.amazonaws.com/ihch_msk_2021.tar.gz.

# Verify data set is 1:1 sample:patient
sample_patient = fread('data_sources/MSK/ihch_msk_2021/data_clinical_sample.txt', skip = 'SAMP')
stopifnot(sample_patient[, unique(.(PATIENT_ID, SAMPLE_ID))][, .N] == uniqueN(sample_patient$PATIENT_ID))

panel_by_sample = fread('data_sources/MSK/ihch_msk_2021/data_gene_panel_matrix.txt')
setnames(panel_by_sample, 'SAMPLE_ID', 'Tumor_Sample_Barcode')

# Split MAF by panel. (In panel_by_sample, the panel name is in the "mutations" field.)
maf = fread("data_sources/MSK/ihch_msk_2021/data_mutations.txt")
maf[panel_by_sample, panel := mutations, on = "Tumor_Sample_Barcode"]

# Two samples are of unknown panel; drop these
maf[is.na(panel), unique(Tumor_Sample_Barcode)]
maf = maf[!is.na(panel)]

by_panel = split(maf, by = "panel", keep.by = F)
output_names = paste0("data_sources/MSK/MSK_2021_IHCH_", names(by_panel), ".maf.gz")
mapply(fwrite, x = by_panel, file = output_names, sep = "\t")


# DZ_EXTENT gives Metastatic disease, Multifocal liver disease, and Solitary liver tumor
# Based on staging manuals, we'll assume cases marked solitary or multifocal are both pM = 0
patient_clinical = fread('data_sources/MSK/ihch_msk_2021/data_clinical_patient.txt', skip = 'PAT')
sample_key = unique(maf[, .(Unique_Patient_Identifier = Tumor_Sample_Barcode)])
sample_key[sample_patient, PATIENT_ID := PATIENT_ID, on = c(Unique_Patient_Identifier = "SAMPLE_ID")]
sample_key[, cca_type := 'IHC']
sample_key[, fluke_status := 'negative'] # per Boerner 2021 Hepatology

# Age is age at diagnosis (per header of data_clinical_patient.txt)
patient_clinical = patient_clinical[, .(age = AGE, sex = SEX, DZ_EXTENT, surv_month = OS_MONTHS, OS_STATUS, PATIENT_ID)]
sample_key = merge.data.table(sample_key, patient_clinical, by = 'PATIENT_ID', all.x = TRUE, all.y = FALSE)
sample_key[, surv_status := fcase(OS_STATUS == '1:DECEASED', 1, OS_STATUS == '0:LIVING', 0)]
sample_key$OS_STATUS = NULL
sample_key[, sex := fcase(sex == 'Female', 'F', sex == 'Male', 'M')]

# Verify DZ_EXTENT values
# table(sample_key$DZ_EXTENT, exclude = NULL)
sample_key[, pM := 0]
sample_key[DZ_EXTENT == 'Metastatic disease', pM := 1]
sample_key$DZ_EXTENT = NULL
setnames(sample_key, 'PATIENT_ID', 'MSK_PATIENT_ID')

# Unique_Patient_Identifier is using MSK sample ID; okay because we have 1:1 samples and patients.
stopifnot(uniqueN(sample_key[, .(MSK_PATIENT_ID, Unique_Patient_Identifier)]) == sample_key[, .N])

# Boerner states "Ninety-two oncogenic fusions were identified in 78/412 patients (19%)."
# This structural variant data file has a roughly similar number of fusion annotations. (But not all
# events that could be fusions are explicitly labeled as fusions in the annotation columns.)
sample_key[, let(fgfr2_fusion = FALSE, other_fusion = FALSE)]
sv = fread("data_sources/MSK/ihch_msk_2021/data_sv.txt")
sv = sv[Sample_Id %in% sample_key$Unique_Patient_Identifier]

# Manual curation: Fusion events are generally described as fusion (or Fusion) in Event_Info.
# Upon reviewing the 29 records that are not described as fusions, three records, all of which contain FGFR,
# are annotated in such a way that they could well be fusions: FGFR2-FLNA, FGFR2-KIAA1217, FGFR2-KRT20. These
# will all be counted as FGFR2 fusions. (One record labeled FGFR2-intragenic will not be counted.)
## sv[! Event_Info %ilike% 'fusion', .(Site1_Hugo_Symbol, Site2_Hugo_Symbol, Event_Info)]
## Note that two passing events gene 1 = gene 2 = FGFR2, but these are separately described as NRAP-FGFR2 fusions.
fgfr2_patients = sv[Site1_Hugo_Symbol == 'FGFR2' | Site2_Hugo_Symbol == 'FGFR2' & 
                  (Event_Info %ilike% 'fusion' | Event_Info %in% c('FGFR2-FLNA', 'FGFR2-KIAA1217', 'FGFR2-KRT20')),
                  unique(Sample_Id)]

# Filter out two events where gene 1 = gene 2 (Comments field describes these as inversions and deletions).
other_patients = sv[Site1_Hugo_Symbol != 'FGFR2' & Site2_Hugo_Symbol != 'FGFR2' & Event_Info %ilike% 'fusion' & 
                      ! Site1_Hugo_Symbol == Site2_Hugo_Symbol, unique(Sample_Id)]

sample_key[fgfr2_patients, fgfr2_fusion := TRUE, on = 'Unique_Patient_Identifier']
sample_key[other_patients, other_fusion := TRUE, on = 'Unique_Patient_Identifier']

# The MSK samples have had silent mutations filtered out (not ideal).
sample_key$excludes_synonymous = TRUE

fwrite(sample_key, 'study_sample_keys/MSK-IMPACT_iCCA_sample_key.txt', sep = "\t")


# Prep gene copy calls
cna = fread('data_sources/MSK/ihch_msk_2021/data_cna.txt')

## According to docs (see
# https://github.com/cBioPortal/cbioportal/blob/v6.0.14/docs/File-Formats.md#discrete-copy-number-data),
# -2 is homozygous deletion, 0 is neutral, 2 is amplification.
stopifnot(all(unlist(cna[, -"Hugo_Symbol"]) %in% c(-2, 0, 2, NA)))
setnames(cna, 'Hugo_Symbol', 'gene')

cna = melt(cna, id.vars = "gene", variable.name = "Unique_Patient_Identifier", 
           value.name = "value", variable.factor = FALSE)
cna[, copy := fcase(value == 0, 'neutral',
                    value == 2, 'amp',
                    value == -2, 'del')]
cna$value = NULL
cna = cna[! is.na(copy)]
cna = cna[Unique_Patient_Identifier %in% sample_key$Unique_Patient_Identifier]

# Deal with genes with outdated or non-standard identifiers.
old_genes = setdiff(cna$gene, ces.refset.hg38$transcripts$gene_name)

# Okay, there is a mix of CDKN2A, CDKN2Ap14ARF, and CDKN2Ap16INK4a.
# Find cases where copy conflicts.
diff_cdkn2a = cna[gene %like% 'CDKN2A', uniqueN(copy), by = 'Unique_Patient_Identifier'][V1 != 1, Unique_Patient_Identifier]

## 6 samples have one isoform called neutral and the other called deleted. We will recode all of
## these as CDKN2A deletions.
cna[diff_cdkn2a, on = "Unique_Patient_Identifier"][gene %like% 'CDKN2A']
# gene Unique_Patient_Identifier value    copy
# <char>                    <char> <int>  <char>
#   1:   CDKN2Ap14ARF           s_WJ_chol_077_T    -2     del
# 2: CDKN2Ap16INK4A           s_WJ_chol_077_T     0 neutral
# 3:   CDKN2Ap14ARF           s_WJ_chol_090_T    -2     del
# 4: CDKN2Ap16INK4A           s_WJ_chol_090_T     0 neutral
# 5:   CDKN2Ap14ARF           s_WJ_chol_108_T    -2     del
# 6: CDKN2Ap16INK4A           s_WJ_chol_108_T     0 neutral
# 7:   CDKN2Ap14ARF           s_WJ_chol_021_T    -2     del
# 8: CDKN2Ap16INK4A           s_WJ_chol_021_T     0 neutral
# 9:   CDKN2Ap14ARF           s_WJ_chol_038_T    -2     del
# 10: CDKN2Ap16INK4A           s_WJ_chol_038_T     0 neutral
# 11:   CDKN2Ap14ARF           s_WJ_chol_049_T     0 neutral
# 12: CDKN2Ap16INK4A           s_WJ_chol_049_T    -2     del

cna[gene %in% c('CDKN2Ap16INK4A', 'CDKN2Ap14ARF'), gene := 'CDKN2A']
cna[Unique_Patient_Identifier %in% diff_cdkn2a & gene == 'CDKN2A', copy := 'del']
cna = unique(cna)

old_genes = old_genes[! old_genes %in% c('CDKN2Ap16INK4A', 'CDKN2Ap14ARF')]

# Plugged remaining genes into syngoportal.org and got some conversions.
conversions = fread('data_sources/MSK/MSK_gene_conversions.txt')
good_conversions = conversions[symbol %in% ces.refset.hg38$transcripts$gene_name]

cna[good_conversions, gene := symbol, on = c(gene = 'query')]

# 11 genes that failed conversion are all universally neutral. They'll be dropped.
still_missing = setdiff(old_genes, good_conversions$query)
stopifnot(cna[gene %in% still_missing, all(copy == 'neutral')])
cna = cna[! gene %in% still_missing]

# We have some new duplicate entries because the converted gene names match some previous records in the CNA calls.
# Therefore, we'll take unique records and remove conflicting calls.
cna = unique(cna)
to_remove =  cna[, .N, by = c('Unique_Patient_Identifier', 'gene')][N > 1, .(Unique_Patient_Identifier, gene)]
cna[to_remove, copy := NA, on = names(to_remove)]
cna = cna[! is.na(copy)]


# We'll subset the other major CNA source (TCGA) to refset coding genes and these three additional
setdiff(cna$gene, ces.refset.hg38$gene_names)
# [1] "RYBP"       "CDKN2B-AS1" "SOX2-OT"   

fwrite(cna, 'gene_copy_calls/MSK_2021_IHCH_gene_copy.txt', sep = "\t")

