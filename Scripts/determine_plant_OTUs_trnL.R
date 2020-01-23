# trnL Plant Sample OTUs
# EKB
# January 22, 2020

# LIBRARIES
library(tidyverse)

# DATA
trnL_reads <- read_csv("Data/SequencedData/Plants/ProcessedData/trnL_reads.csv")
trnL_totals <- read_csv("Data/SequencedData/Plants/ProcessedData/trnL_totals.csv")
sample_id <- read_csv("Data/CollectionData/vial_id.csv") %>% 
  filter(sample_type == "plant")

# Get only plant vouchers
trnL_reads <- filter(trnL_reads, SampleID %in% sample_id$vial_id)
trnL_totals <- filter(trnL_totals, SampleID %in% sample_id$vial_id)

# Make Histograms
