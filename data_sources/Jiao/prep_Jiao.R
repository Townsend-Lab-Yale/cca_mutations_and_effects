library(readxl)
library(data.table)
library(cancereffectsizeR)
library(ces.refset.hg38)

# Read in the Jiao calls from the Excel file. They use hg18 coordinates.
maf = as.data.table(read_excel('data_sources/Jiao/original/Jiao_S4.xlsx', skip = 2))
maf = maf[1:(.N - 1)] # last row is footer material

# Parse chr/start/ref/alt from Jiao's format (e.g., g.chr9:106660684C>A)
maf[, c('Chromosome', 'change') := tstrsplit(`Nucleotide (genomic)`, split = ':')]
maf[, Chromosome := sub('^g\\.chr', '', Chromosome)]
maf[, Start_Position := as.numeric(sub('(^\\d+).*', '\\1', change))]
maf[, ref_alt := sub('.*\\d', '', change)]
maf[, c('Reference_Allele', 'Tumor_Allele') := tstrsplit(ref_alt, split = '>')]

# Substitutions (including dinculeotide substitutions) have been handled. Indels need their ref/alt fixed.
# For indels the above manipluations lead to NA Tumor_Allele and everything in the ref field (e.g., delG or insGGC)
maf[`Mutation type*` == 'DEL', let(Reference_Allele = sub('del', '', Reference_Allele),
                                   Tumor_Allele = '-')]
maf[`Mutation type*` == 'INS', let(Reference_Allele = '-', Tumor_Allele = sub('ins', '', Reference_Allele))]

maf[, Tumor_Sample_Barcode := paste0('jiao_', `Tumor Sample`)]

# Lift to hg38. A handful of records will get dropping for failing liftOver or no longer having correct reference allele.
lifted = preload_maf(maf = maf, refset = 'ces.refset.hg38', chain_file = 'reference/chains/hg18ToHg38.over.chain')

# ~0.35% of records had liftOver problems.
stopifnot(lifted[, mean(! is.na(problem) & problem %in% c('failed_liftOver', 'reference_mismatch',
                                                                      'not_variant', 'duplicate_record_after_liftOver'))] - .0035 < 1e-4)
lifted = lifted[is.na(problem)]


# Remove the GBC sample (discernable from sample identifiers).
lifted = lifted[Unique_Patient_Identifier %like% 'CHOL']
fwrite(lifted, 'data_sources/Jiao/jiao.maf.gz', sep = "\t")

sample_key = data.table(Unique_Patient_Identifier = unique(lifted$Unique_Patient_Identifier), cca_type = 'IHC')
sample_key = sample_key[order(Unique_Patient_Identifier)]
fwrite(sample_key, 'study_sample_keys/jiao_sample_key.txt', sep = "\t")





  