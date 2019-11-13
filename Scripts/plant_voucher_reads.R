# Get Plant Voucher Reads
# 12 November 2019
# EKB

# LIBRARIES #
library(tidyverse)

# FILES #



trnL_fall17 <- read_csv("Data/SequencedData/Plants/RawData/trnL_Fall2017.csv")

# DATA PREP #

# Plant Vial IDs
vial_id <- read_csv("Data/CollectionData/vial_id.csv")
plant_vial_id <- vial_id %>% filter(sample_type == "plant")

# trnL Fall 2016
trnL_fall16 <- read_csv("Data/SequencedData/Plants/RawData/trnL_Fall2016.csv")

# trnL Spring 2017
trnL_spring17 <- read_csv("Data/SequencedData/Plants/RawData/trnL_Spring2017.csv")