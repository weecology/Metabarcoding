# trnL data (read in from plots_for_manuscript.r)

reads <- read_csv("Data/SequencedData/Plants/ProcessedData/trnL_reads.csv")
totals <- read_csv("Data/SequencedData/Plants/ProcessedData/trnL_totals.csv")
samples <- read_csv("Data/CollectionData/fecal_sample_collection.csv")

trnL_all <- filter_reads_data_trnL(samples,
                       reads,
                       totals,
                       reads_min = 2000,
                       rel_reads_min = 0.01)

trnL_all <- left_join(trnL_all[[3]], 
                      trnL_all[[1]] %>% select(vial_barcode, group), 
                      by = c("SampleID" = "vial_barcode"))
trnL_all <- trnL_all %>% 
  group_by(group) %>% 
  filter(group != 'CP: KR Exclosure') %>% 
  distinct_at(.vars = "OTU")

CP_C_OTUs <- filter(trnL_all, group == 'CP: Control') %>% 
  ungroup() %>% 
  select(OTU)
Dipo_OTUs <- filter(trnL_all, group == 'K-Rat') %>% 
  ungroup() %>% 
  select(OTU)

overlap_trnL <- as.data.frame(count(intersect(CP_C_OTUs, Dipo_OTUs)))
total_trnL <- length(unique(trnL_all$OTU))

percent_overlap_trnL <- overlap_trnL/total_trnL
percent_overlap_trnL

(percent_Dipo_diet <- overlap_trnL/as.data.frame(count(Dipo_OTUs)))

# ITS2 data 

reads <- read_csv("Data/SequencedData/Plants/ProcessedData/ITS2_reads_WeeTU.csv")
totals <- read_csv("Data/SequencedData/Plants/ProcessedData/ITS2_totals.csv")
samples <- read_csv("Data/CollectionData/fecal_sample_collection.csv")

# get plant OTUs only
reads <- filter(reads, WTU.clade1 == 4) %>% 
  select(OTU:DataFrame)

ITS2_all <- filter_reads_data_ITS2(samples,
                                   reads,
                                   totals,
                                   reads_min = 2000,
                                   rel_reads_min = 0.01)

ITS2_all <- left_join(ITS2_all[[3]], 
                      ITS2_all[[1]] %>% select(vial_barcode, group), 
                      by = c("SampleID" = "vial_barcode"))
ITS2_all <- ITS2_all %>% 
  group_by(group) %>% 
  filter(group != 'CP: KR Exclosure') %>% 
  distinct_at(.vars = "OTU")

CP_C_OTUs <- filter(ITS2_all, group == 'CP: Control') %>% 
  ungroup() %>% 
  select(OTU)
Dipo_OTUs <- filter(ITS2_all, group == 'K-Rat') %>% 
  ungroup() %>% 
  select(OTU)

overlap_ITS2 <- as.data.frame(count(intersect(CP_C_OTUs, Dipo_OTUs)))
total_ITS2 <- length(unique(ITS2_all$OTU))

percent_overlap_ITS2 <- overlap_ITS2/total_ITS2
percent_overlap_ITS2

(percent_Dipo_diet <- overlap_ITS2/as.data.frame(count(Dipo_OTUs)))
