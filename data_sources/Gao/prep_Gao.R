library(data.table)
meta = fread('data_sources/Gao/gao_sequencing_data_accessions.txt')
meta = meta[sample_type == 'tumor'][, .(Sample_Name, age, sex, disease)]
meta[, sex := fcase(sex == 'female', 'F', sex == 'male', 'M')]
meta[, Unique_Patient_Identifier := paste0('gao_', Sample_Name)]
stopifnot(all(meta$disease == "intrahepatic cholangiocarcinoma"))
meta[, cca_type := 'IHC']
meta[, pM := 0] # Gao supplementary table 1

# Paper says, "These well-characterized patients received no prior anticancer treatments and were
# confirmed to have no background liver or biliary diseases."
meta[, fluke_status := 'negative']
meta = meta[, .(Unique_Patient_Identifier, cca_type, fluke_status, age, sex, pM)]
fwrite(meta, 'study_sample_keys/gao_sample_key.txt', sep = "\t")
