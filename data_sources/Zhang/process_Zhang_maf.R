library(data.table)
library(cancereffectsizeR)

# MAF file acquired directly from authors
zhang_maf_file = 'data_sources/Zhang/original/Zhang_from_authors.maf.gz'
zhang_maf = fread(zhang_maf_file, quote = "")

sample_key = unique(zhang_maf[, .(Unique_Patient_Identifier = paste0('zhang_', Tumor_Sample_Barcode), 
                                  cca_type = type)])
sample_key[cca_type == 'IHCC', cca_type := 'IHC']
sample_key[cca_type == 'Perihilar', cca_type := 'PHC']
stopifnot(all(sample_key$cca_type %in% c('PHC', 'IHC')))
fwrite(sample_key, 'study_sample_keys/zhang_sample_key.txt', sep = "\t")

# fix a column name and restrict to standard contigs
setnames(zhang_maf, 'Start_position', 'Start_Position')
zhang_maf = zhang_maf[Chromosome %in% c(1:22, 'X', 'Y')]
zhang_maf[, Tumor_Sample_Barcode := paste0('zhang_', Tumor_Sample_Barcode)]
zhang_maf = zhang_maf[, .(Tumor_Sample_Barcode, Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2)]
fwrite(zhang_maf, 'data_sources/Zhang/zhang.maf.gz', sep = "\t")



