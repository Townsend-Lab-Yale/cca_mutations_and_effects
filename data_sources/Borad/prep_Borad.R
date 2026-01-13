library(data.table)
library(ces.refset.hg19)

borad = fread("data_sources/Borad/Borad_PLOS_Genetics_S1_conv.txt")
borad = borad[1:(.N - 1)] # last row has footer content
borad = borad[, .(Chromosome = Chr, Start_Position = hg19_position, 
                  Reference_Allele = Reference_allele, Tumor_Allele = Alternate_allele,
                  Tumor_Sample_Barcode = paste0('borad_', Patient))]

nt = c('A', 'C', 'G', 'T')
genome = ces.refset.hg19$genome

# confirm correct reference alleles for deletions
stopifnot(borad[! Reference_Allele %in% nt, 
                identical(as.character(getSeq(genome, GRanges(paste0(Chromosome, ':', Start_Position, '-', 
                                                                     Start_Position + nchar(Reference_Allele) - 1)))),
                          Reference_Allele)])

# Ref/Alt alleles are VCF-style (as in, GGGA>G rather than GGA>-), so we
# need to increment start position and fix alleles
borad[! Reference_Allele %in% nt, 
      c("Start_Position", "Reference_Allele", "Tumor_Allele") := 
          .(Start_Position + 1, substr(Reference_Allele, 2, nchar(Reference_Allele)), '-')]

# Fix one indel record that is in a unique format
borad[Tumor_Allele == '<DEL>', Tumor_Allele := '-']

# For inserts, Start_Position stays the same, but change style (T/TG becomes -/G)
borad[! Tumor_Allele %in% c('-', nt), c("Reference_Allele", "Tumor_Allele") :=
        .('-', substr(Tumor_Allele, 2, nchar(Tumor_Allele)))]


# Split by exome capture kit as given in Borad methods.
v1 = borad[Tumor_Sample_Barcode %in% c('borad_1', 'borad_3')]
nextera = borad[Tumor_Sample_Barcode == 'borad_2']
truseq = borad[Tumor_Sample_Barcode %in% c('borad_4', 'borad_5')]
v4_utr = borad[Tumor_Sample_Barcode == 'borad_6']

fwrite(v1, 'data_sources/Borad/borad_SureSelectV1.maf.gz', sep = "\t")
fwrite(nextera, 'data_sources/Borad/borad_Nextera.maf.gz', sep = "\t")
fwrite(truseq, 'data_sources/Borad/borad_Truseq.maf.gz', sep = "\t")
fwrite(v4_utr, 'data_sources/Borad/borad_SureSelectV4-UTR.maf.gz', sep = "\t")

# Note: Borad sample key file (see study_sample_keys) produced manually based on information in the Borad paper.
patient_data = transpose(as.data.table(readxl::read_excel('data_sources/Borad/borad_patient_data_from_Fig5.xlsx', skip = 2)),
                         keep.names = 'Unique_Patient_Identifier', make.names = 1)
patient_data[, Unique_Patient_Identifier := sub('Patient ', 'borad_', Unique_Patient_Identifier)]
patient_data[`Location of Primary Tumor` %like% 'Intrahepatic', cca_type := 'IHC']
stopifnot(all(patient_data$`Liver fluke` == 'No'))
patient_data[, fluke_status := 'negative']

# Patient 2 is listed as Stage IV, probably under AJCC edition 7. However, the patient has only regional (abdominal) lymph
# node metastasis, which is not distant metastasis (and wouldn't be stage IV under 8th edition).
patient_data[, pM := ifelse(Stage == 'IV', 1, 0)]
patient_data[Unique_Patient_Identifier == 'borad_2', pM := 0]
patient_data[, let(age = `Age (years)`,
                   sex = Gender)]

# Recode survival data; drop plus signs on some month survival
patient_data[, surv_status := ifelse(`Survival Status` == 'Alive', 0, 1)]
patient_data[, surv_month := as.numeric(sub('\\+$', '', `Survival Duration from biopsy (months)`))]


# Borad's paper's Table 7 reports FGFR2 fusions in borad_4, borad_5, borad_6.
patient_data[, fgfr2_fusion := Unique_Patient_Identifier %in% c('borad_4', 'borad_5', 'borad_6')]


borad_key = patient_data[, .(Unique_Patient_Identifier, cca_type, age, sex, pM, surv_status, surv_month, 
                             fluke_status, fgfr2_fusion)]

fwrite(borad_key, 'study_sample_keys/borad_sample_key.txt', sep = "\t")
