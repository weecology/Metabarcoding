# Prep trnL data from Jonah Inc.
# Ellen Bledsoe
# Jan 2017

########################
# LIBRARIES

library(dplyr)
<<<<<<< HEAD
library(stringr)

########################
# LOAD FILES

# load the files you need
trnL <- read.csv("./", header = TRUE)
samples <- read.csv("./", header = TRUE)
=======

# set working directory
setwd("~bleds22e/Documents/Git/Metagenomics/")
>>>>>>> database_creation

########################
# LOAD FILES

trnL_fecal <- read.csv("./Plants/trnL_fecal_samples.csv", header = TRUE)
fecal <- read.csv("./CollectionData/fecal_sample_collection.csv", header = TRUE)
trnL_voucher <- read.csv("./Plants/trnL_voucher_specimens.csv", header = TRUE)
voucher <- read.csv("./CollectionData/plant_voucher_collection.csv", header = TRUE)
trnL_link_file <- read.csv("./Plants/trnL_link_file.csv", header = TRUE)

#=====================================================================================
### FULL trnL LINK FILE ###

# rename columns
colnames(trnL_link_file) <- c("OTU_cluster", "Type", "Taxa")

for(this_level in c('d','k','p','c','o','f','g','s')){
  # separate taxa into columns
  step_one = sapply(strsplit(as.character(trnL_link_file$Taxa), paste0(this_level,':')), '[', 2)
  step_two = sapply(strsplit(step_one, ','), '[', 1)
  trnL_link_file[ ,this_level] = step_two
}

# remove trailing semicolon from species column
trnL_link_file$s = sapply(strsplit(trnL_link_file$s, ';'), '[', 1)

# rename columns
trnL_link_file <- select(trnL_link_file, OTU_cluster, Type, d, k, p, c, o, f, g, s) %>% 
  rename(Family = o, Genus = g, SciName = s)

# remove '*' from Type column
trnL_link_file$Type <- str_replace(trnL_link_file$Type, '[*]', NA_character_)

# write CSV file
# write.csv(trnL_link_file, "./SequencedData/Plants/trnL_taxa_link_file.csv", row.names = FALSE)

#=====================================================================================
### FECAL DATA ###

# remove unnecessary columns summary rows
trnL_fecal <- trnL_fecal[-(c(452:453)), c(1,8:54)]

# turn into flat file, keeping only non-zero numbers
trnL_fecal <- tidyr::gather(trnL_fecal, "Sample", "Proportion", 2:48) %>% 
  filter(Proportion != 0)

# write CSV file
# write.csv(trnL_fecal, "./SequencedData/Plants/trnL_fecal_data.csv", row.names = FALSE)

#=====================================================================================
### VOUCHER DATA ###

# remove unnecessary columns summary rows
trnL_voucher <- trnL_voucher[-c(14205:14206), -c(2,108)]

# turn into flat file with only non-zero values
trnL_voucher <- tidyr::gather(trnL_voucher, "Sample", "Reads", 2:106) %>% 
  filter(Reads != 0) 

# split VialBarcode column and keep only vial number
trnL_voucher <- tidyr::separate(trnL_voucher, Sample, c("Sample", "Wisely"), remove = TRUE) %>%
  select(-Wisely)

# write CSV file
# write.csv(trnL_voucher, "./SequencedData/Plants/trnL_voucher_data.csv", row.names = FALSE)
