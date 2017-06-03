# Prep trnL data from Jonah Inc.
# Ellen Bledsoe
# Jan 2017

########################
# LIBRARIES

library(dplyr)
library(stringr)

########################
# LOAD FILES

# set your working directory
setwd("~bleds22e/Documents/Git/Metagenomics/")

# load the files you need
trnL <- read.csv("./SequencedData/Trap_Bait_Test/trnL_trap_and_bait_test.csv", header = TRUE)
fecal <- read.csv("./CollectionData/fecal_sample_collection.csv", header = TRUE)
plants <- read.csv("./CollectionData/plant_voucher_collection.csv", header = TRUE)

trnL_voucher_data <- read.csv("./SequencedData/Plants/trnL_voucher_data.csv", header = T)

########################
# CLEAN DATA

# reads dataframe
reads <- trnL[,-c(2,100)]
reads <- tidyr::gather(reads, "Sample", "Reads", 2:98) %>% 
  filter(Reads != 0)
reads <- tidyr::separate(reads, Sample, into = c("Sample", "Wisely")) %>%
  select(-Wisely)

# vouchers dataframe
vouchers <- select(plants, vial_barcode, year) %>% 
  rename(Sample = vial_barcode) %>% 
  select(-year)
vouchers <- semi_join(reads, vouchers, by = "Sample")

new_file <- rbind(trnL_voucher_data, vouchers)
#write.csv(new_file, "SequencedData/Plants/trnL_voucher_data.csv", row.names = FALSE)

# fecal sample dataframe

fecal <- anti_join(all_ITS, vouchers, by = "Sample")

#write.csv(fecal, "SequencedData/Plants/ITS_fecal_data.csv", row.names = F)