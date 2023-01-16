# preprocess data
# 
# ### rebuild dds
message("building dds for EIS dataset")
dds.EIS <- ddsFromFeatureCounts(file.path(CNT_FILES$path, CNT_FILES$EIS.cnts), 
                                file.path(CNT_FILES$path, CNT_FILES$EIS.conditions), 
                                source = "tximport", 
                                gffa3 = file.path(A_FILES$path, A_FILES$gff3), 
                                design = ~ Genotype + Treatment + Genotype:Treatment)

# message("building dds for EIS")
# dds.EIS <- ddsFromFeatureCounts(file.path(CNT_FILES$path, CNT_FILES$EISd75.cnts), 
#                                    file.path(CNT_FILES$path, CNT_FILES$EISd75.conditions),
#                                    source = "tximport", 
#                                    gene_info = file.path(A_FILES$path, A_FILES$annotation), 
#                                    design = ~ Genotype + Treatment + Genotype:Treatment)

message("building dds for LCM dataset")
dds.LCM <- ddsFromFeatureCounts(file.path(CNT_FILES$path, CNT_FILES$LCM.cnts), 
                                file.path(CNT_FILES$path, CNT_FILES$LCM.conditions),
                                source = "tximport", 
                                gffa3 = file.path(A_FILES$path, A_FILES$gff3), 
                                design = ~ Genotype + Treatment + Genotype:Treatment)


EIS.RES <- importResults(EIS_RES_FILES)
LCM.RES <- importResults(LCM_RES_FILES)

GENE_FAMILY <- read_delim(file.path(A_FILES$path, A_FILES$genefamily), 
                          delim = "\t", quote = "", col_types = cols()) %>%
              select(Family = Gene_Family, everything())%>%
              filter(!Genomic_Locus_Tag == "NULL") %>%
              mutate(Sub_Family = ifelse(Sub_Family=="NULL", Family, Sub_Family))

message("done...")