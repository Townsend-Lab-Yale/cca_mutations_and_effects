library(data.table)
library(readxl)
library(BSgenome.Hsapiens.UCSC.hg19)


genome = getBSgenome(BSgenome.Hsapiens.UCSC.hg19)
seqlevelsStyle(genome) = 'NCBI' # warning okay

# Read in file containing Zou supplementary data 3
s3 = as.data.table(read_excel('data_sources/Zou/original/41467_2014_BFncomms6696_MOESM1279_ESM.xlsx', skip = 1))
s3 = s3[2:.N] # first row is a remnant from double-row header
s3 = s3[, .(Tumor_Sample_Barcode = paste0('zou_', Patient), Type,
            Chromosome = Chrom, Start_Position = Pos, Reference_Allele = Ref, Tumor_Seq_Allele2 = Var)]

# remove mitochondrial
s3 = s3[Chromosome != 'MT']

# Deletions appear to have Ref/Alt like -/-AAG when they should be AAG/- for MAF format.
# Verify that this is correct by checking reference alleles at sites, then change format.
del_check = s3[Type == 'DEL', .(seqnames = Chromosome, start = Start_Position, 
                                end = Start_Position + nchar(Tumor_Seq_Allele2) - 2, 
                                assumed_ref = substr(Tumor_Seq_Allele2, 2, nchar(Tumor_Seq_Allele2)))]
del_alleles = as.character(getSeq(genome, makeGRangesFromDataFrame(del_check)))
stopifnot(identical(del_check$assumed_ref, del_alleles))

# The check passed, so fix the format accordingly
s3[Type == 'DEL', c("Reference_Allele", "Tumor_Seq_Allele2") := .(substr(Tumor_Seq_Allele2, 2, nchar(Tumor_Seq_Allele2)),
                                                                  '-')]

# Insertions are like C/+GT instead of -/GT. (Always one reference base, and a plus sign in the ALT allele.)
s3[Type == 'INS', c("Reference_Allele", "Tumor_Seq_Allele2") := .('-', 
                                                                  substr(Tumor_Seq_Allele2, 2, nchar(Tumor_Seq_Allele2)))]

s3 = s3[, -"Type"]
fwrite(s3, 'data_sources/Zou/zou.maf.gz', sep = "\t")

# Zou paper says "patients studied here do not have liver fluke infection."
sample_key = data.table(Unique_Patient_Identifier = unique(s3$Tumor_Sample_Barcode), cca_type = 'IHC', 
                        fluke_status = 'negative')

# Zou patient info is in a supplemental file that uses different identifiers than the MAF!
# Luckily, there is ultimately a 1:1 correspondence between the numerical parts of the identifiers.
# For some reason, a single sample is on its own separate sheet.
sheet1 = as.data.table(readxl::read_excel('data_sources/Zou/original/Zou_supp_data_1.xls', skip = 1))
sheet2 = as.data.table(readxl::read_excel('data_sources/Zou/original/Zou_supp_data_1.xls', sheet = 2))

# match up other sheet
setnames(sheet2, c('survival time (m)', 'Age', 'Sex'), c('survival(m)', 'age', 'sex'))

# Fix column types to match so sheets can go to rbind.
sheet1$`following time` = as.character(sheet1$`following time`)
sheet2$`time of operation` = as.character(sheet2$`time of operation`)
sheet2$`time of recurrence` = as.character(sheet2$`time of recurrence`)
sheet2$`death time` = as.character(sheet2$`death time`)
sheet2$`following time` = as.character(sheet2$`following time`)

patient_info = rbind(sheet1, sheet2, fill = T)


# All but four samples have an ID value containing CT followed by the sample number seen in the MAF
setnames(patient_info, 'Original Sample Name...2', 'tmp_id')
patient_info[tmp_id %like% 'CT\\d+', new_id := paste0('zou_', sub('.*CT(\\d+).*', '\\1', tmp_id))]

# Remaining four samples have their sample numbers (which are the only sample numbers that also contain hyphens)
# encoded in a different format.
patient_info[! tmp_id %like% 'CT' & tmp_id %like% '\\(8', new_id := paste0('zou_', sub('.*\\(8-(.*)\\)$', '8_\\1', tmp_id))]

# new_id is uniquely identifying
stopifnot(uniqueN(patient_info$new_id) == patient_info[, .N])

# We have recovered sample info for 102 of 103 Zou samples.
stopifnot(length(intersect(sample_key$Unique_Patient_Identifier, patient_info$new_id)) == 102,
          sample_key[, .N] == 103, patient_info[, .N] == 103)

# Confirmed that survival is coded as 1 = deceased, 0 = alive; and extrahepatic metastasis also binary.
patient_info = patient_info[, .(Unique_Patient_Identifier = new_id, sex, age, surv_month = `survival(m)`, 
                                surv_status = `survival status`, pM = `extrahepatic metastasis`)]

# According to the key on sheet 3 of the Excel file
patient_info[, sex := fcase(sex == 0, 'F', sex == 1, 'M')]

sample_key = merge.data.table(sample_key, patient_info, all.x = TRUE, all.y = FALSE, 
                              by = 'Unique_Patient_Identifier')
fwrite(sample_key, 'study_sample_keys/zou_sample_key.txt', sep = "\t")
