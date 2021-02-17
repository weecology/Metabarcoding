# read FNA files
# EKB
# 2/9/2017

# to run this script, you will need to download bioconductor first:
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")

BiocManager::install("Biostrings")
BiocManager::install("annotate")

#===== LIBRARIES =====#

library(Biostrings)

#===== READ IN FILES =====#

trnL_fall2016 <- readDNAStringSet("../../Portal/PORTAL_primary_data/DNA/Results_Jonah/Plants/Fall2016/trnL_refseqs_022417.fna")
its_fall2016 <- readDNAStringSet("../../Portal/PORTAL_primary_data/DNA/Results_Jonah/Plants/Fall2016/ITS2_merged.prtrim.upf.filt.derep.mc2.repset.fna")

unkn_fall2017 <- readDNAStringSet("../../Portal/PORTAL_primary_data/DNA/Results_Jonah/Plants/Fall2017/R1.170.prtrim.upf.filt.derep.mc2.repset.fna")
ref_samples <- readDNAStringSet("../../Portal/PORTAL_primary_data/DNA/Results_Jonah/Plants/reference_samples/ref_seq_otus.fna")

spring2017 <- readDNAStringSet("../../Portal/PORTAL_primary_data/DNA/Results_Jonah/Plants/Spring2017/r1.prtrim.upf.filt.derep.mc2.repset.fna")
corrected_OTUs_spring2017 <- readDNAStringSet("../../Portal/PORTAL_primary_data/DNA/Results_Jonah/Plants/Spring2017/correctedOTUs_r1.prtrim.upf.filt.derep.mc2.repset.fna")

# ===== MAKE DATAFRAMES =====#
OTU_its <- names(its)
sequence_its <- paste(its)
df_its <- data.frame(OTU_its, sequence_its)

OTU_trnL <- names(trnL)
sequence_trnL <- paste(trnL)
df_trnL <- data.frame(OTU_trnL, sequence_trnL)


#===== WRITE CSV FILES =====#

write.csv(df_its, "data/its_otus.csv")
write.csv(df_trnL, "Plants/trnL_from_fna.csv", row.names=FALSE)




