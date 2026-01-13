maf = fread('data_sources/Chan-on//Chan-on_CCA.maf.gz')

# From Chan-on supplementary information
v1_samples = c("chan-on_72600277T", "chan-on_33587479T", "chan-on_21914187T", "chan-on_02000861T", 
               "chan-on_29150572T", "chan-on_17231203T", "chan-on_0715T", "chan-on_0824T", 
               "chan-on_0930T", "chan-on_1202T")

v4_utr_samples =  c("chan-on_00990196T", "chan-on_02000123T", "chan-on_23474504T", 
                "chan-on_77071507T", "chan-on_3001T")

stopifnot(all(unique(maf$Tumor_Sample_Barcode) %in% c(v1_samples, v4_utr_samples)))

v1 = maf[Tumor_Sample_Barcode %in% v1_samples]
v4_utr = maf[Tumor_Sample_Barcode %in% v4_utr_samples]

fwrite(v1, 'data_sources/Chan-on/chan-on_SureSelectV1.maf.gz', sep = "\t")
fwrite(v4_utr, 'data_sources/Chan-on/chan-on_SureSelectV4-UTR.maf.gz', sep = "\t")
