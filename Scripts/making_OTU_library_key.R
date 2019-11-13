# Making OTU Library Key
# 12 November 2019
# EKB

# This script is for making an OTU library for our plant vouchers. 
#   - there will be separate files for trnL and ITS2
#   - each file will contain:
#       - the OTU from Jonah Ventures
#       - the family, genus, and species of the OTU
#       - the WeeTU for filtering by family, genus, or species

# Current questions:
#   - Are we basing this only on plant vouchers that were sent in?
#   - OR from Portal plant list?
#   - OR are we going through OTUs from Jonah Ventures that make family, genus,
#     species that we have at the site?
#   - OR the OTUs that have number/proportion of reads above a certain threshold?

# TBD before really digging in: do I have data from all the plant vouchers?
#   - pull out voucher data first? Then convert to reads or proportions or not
#     even bother with that part?
#   - might want to source some code for pulling out sample numbers, so you can 
#     use new vial_id csv for pulling out plant voucher samples

#==============================================================================

# LIBRARIES #
library(tidyverse)

# PREP MATERIALS #

# Prep trnL reference csv into functional dataframe #

trnL <- read_csv("Data/SequencedData/Plants/RawData/JV13trnlClosedRef.csv")
colnames(trnL) <- c("OTU_cluster", "Seed", "ConsensusLineage")

# split ConsesusLineage column into usable columns
for(this_level in c('d','k','p','c','o','f','g','s')){
  # separate taxa into columns
  step_one <- sapply(strsplit(as.character(trnL$ConsensusLineage), paste0(this_level,':')), '[', 2)
  step_two <- sub(";", "", step_one)
  step_three <- sapply(strsplit(step_two, ','), '[', 1)
  trnL[,this_level] <- step_three
}

# split species column into two columns and filter for SEED rows
trnL <- trnL %>%
  select(., -g) %>% 
  tidyr::separate(s, c("g", "s"), " ") %>% 
  filter(Seed == "SEED")