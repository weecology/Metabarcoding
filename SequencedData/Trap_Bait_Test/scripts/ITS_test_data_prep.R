# Prep trap & bait test ITS2 data from Jonah Inc.
# Ellen Bledsoe
# June 2017

########################
# LIBRARIES

library(dplyr)

########################
# LOAD FILES

# set working directory
setwd("~bleds22e/Documents/Git/Metagenomics/")

# load files
ITS <- read.csv("./SequencedData/Trap_Bait_Test/original_data/ITS_blasted_trap_and_bait_test.csv", header = TRUE, na.strings = "")
fecal <- read.csv("./CollectionData/fecal_sample_collection.csv", header = TRUE)
plants <- read.csv("./CollectionData/plant_voucher_collection.csv", header = TRUE)

########################
# CLEAN DATA

# remove Sum column
ITS <- select(ITS, -Sum)

# Pull out ConsensusLineage

its_test_taxa <- select(ITS, OTU.ID, ConsensusLineage)

for(this_level in c('k','p','c','o','f','g','s')){
  # separate taxa into columns
  step_one=sapply(strsplit(as.character(its_test_taxa$ConsensusLineage), paste0(this_level,'__')), '[', 2)
  step_two=sapply(strsplit(step_one, ';'), '[', 1)
  its_test_taxa[,this_level]=step_two
}

# make taxa link file
its_test_taxa <- select(its_test_taxa, OTU.ID, k:s) %>% 
  rename(Family = f, Genus = g, SciName = s)
write.csv(its_test_taxa, "SequencedData/Trap_Bait_Test/ITS_test_taxa_link_file.csv", row.names = FALSE)

# Restructure

ITS <- select(ITS, -ConsensusLineage)

ITS_reads <- tidyr::gather(ITS, "Sample", "Reads", 2:98) %>%
  filter(Reads != 0)
ITS_reads <- tidyr::separate(ITS_reads, Sample, into = c("Sample", "Wisely")) %>%
  select(-Wisely)

# vouchers dataframe
vouchers <- select(plants, vial_barcode, year) %>% 
  rename(Sample = vial_barcode) %>% 
  select(-year)
vouchers_its <- semi_join(ITS_reads, vouchers, by = "Sample")
#write.csv(vouchers_its, "SequencedData/Trap_Bait_Test/test_ITS_voucher_data.csv", row.names = FALSE)

# fecal sample dataframe
fecal <- anti_join(ITS_reads, vouchers, by = "Sample")
#write.csv(fecal, "SequencedData/Trap_Bait_Test/test_ITS_fecal_data.csv", row.names = F)
