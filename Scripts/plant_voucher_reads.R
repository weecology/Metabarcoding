# Get Plant Voucher Reads
# 12 November 2019
# EKB

# LIBRARIES #
library(tidyverse)

#==============================================================================#
# DATA PREP: trnL #

## Plant Vial IDs ##
vial_id <- read_csv("Data/CollectionData/vial_id.csv")
plant_vial_id <- vial_id %>% filter(sample_type == "plant")
fecal_vial_id <- vial_id %>% filter(sample_type == "fecal")

## trnL Fall 2016 ##

trnL_fall16 <- read_csv("Data/SequencedData/Plants/RawData/trnL_Fall2016.csv")

# remove extra columns
trnL_fall16 <- trnL_fall16[, c(1, 8:55)]
colnames(trnL_fall16)[1] <- "OTU"

# remove Sum column and total reads row
trnL_fall16_total_reads <- trnL_fall16[nrow(trnL_fall16),]
trnL_fall16_OTUsum <- trnL_fall16[, c(1, length(trnL_fall16))]
trnL_fall16 <- trnL_fall16[-nrow(trnL_fall16), -length(trnL_fall16)]

## trnL Spring 2017 ##

trnL_spring17 <- read_csv("Data/SequencedData/Plants/RawData/trnL_Spring2017.csv")

# fix naming (i.e., remove Wisely from vial IDs)
colnames(trnL_spring17) <- sub(".Wisely", "", colnames(trnL_spring17))
colnames(trnL_spring17)[1] <- "OTU"

# remove extra columns and rows
trnL_spring17 <- trnL_spring17[,-2]
trnL_spring17$Sum <- rowSums(trnL_spring17[,-1])
trnL_spring17 <- trnL_spring17[!(trnL_spring17$Sum == 0),]
trnL_spring17_OTUsum <- trnL_spring17[, c(1, length(trnL_spring17))]
trnL_spring17_total_reads <- trnL_spring17[nrow(trnL_spring17),]
trnL_spring17 <- trnL_spring17[-nrow(trnL_spring17),-length(trnL_spring17)]

## trnL Fall 2017 ##

trnL_fall17 <- read_csv("Data/SequencedData/Plants/RawData/trnL_Fall2017.csv")

# fix naming (i.e., remove Bledsoe from vial IDs)
colnames(trnL_fall17) <- sub(".Bledsoe", "", colnames(trnL_fall17))
colnames(trnL_fall17)[1] <- "OTU"

# remove extra columns and rows
trnL_fall17 <- trnL_fall17[-(nrow(trnL_fall17)-1),-2]
trnL_fall17 <- trnL_fall17[!(trnL_fall17$Sum == 0),]
trnL_fall17_total_reads <- trnL_fall17[nrow(trnL_fall17),]
trnL_fall17_OTUsum <- trnL_fall17[, c(1, length(trnL_fall17))]
trnL_fall17 <- trnL_fall17[-(nrow(trnL_fall17)), -length(trnL_fall17)]

## trnL Vouchers ##

vouchers <- read_csv("Data/SequencedData/Plants/RawData/trnL_refsamples.csv")

# remove extra columns
vouchers <- vouchers[, -2]
colnames(vouchers)[1] <- "OTU"
colnames(vouchers) <- sub(".Wisely", "", colnames(vouchers))

vouchers$Sum <- rowSums(vouchers[,-1])
vouchers_OTUsum <- vouchers[, c(1, length(vouchers))]
vouchers_total_reads <- vouchers[nrow(vouchers),]
vouchers <- vouchers[-nrow(vouchers),-length(vouchers)]

#==============================================================================#
# COMBINE DATASETS: trnL #

# convert spring 2017 to proportions
trnL_spring17_OTUs <- trnL_spring17[,1]
trnL_spring17 <- trnL_spring17 %>% 
  column_to_rownames(var = "OTU")
trnL_spring17_total_reads <- trnL_spring17_total_reads[,-c(1, length(trnL_spring17_total_reads))]

trnL_spring17 <- mapply('/', trnL_spring17, trnL_spring17_total_reads)
trnL_spring17 <- as_tibble(trnL_spring17)
trnL_spring17 <- cbind(trnL_spring17_OTUs, trnL_spring17)

# convert fall 2017 to proportions
trnL_fall17_OTUs <- trnL_fall17[,1]
trnL_fall17 <- trnL_fall17 %>% 
  column_to_rownames(var = "OTU")
trnL_fall17_total_reads <- trnL_fall17_total_reads[, -c(1, length(trnL_fall17))]

trnL_fall17 <- mapply('/', trnL_fall17, trnL_fall17_total_reads)
trnL_fall17 <- as_tibble(trnL_fall17)
trnL_fall17 <- cbind(trnL_fall17_OTUs, trnL_fall17)

# convert vouchers to proportions
vouchers_OTUs <- vouchers[,1]
vouchers <- vouchers %>% 
  column_to_rownames(var = "OTU")
vouchers_total_reads <- vouchers_total_reads[,-c(1, length(vouchers_total_reads))]

vouchers <- mapply('/', vouchers, vouchers_total_reads)
vouchers <- as_tibble(vouchers)
vouchers <- cbind(vouchers_OTUs, vouchers)

# combine dataframes
trnL_props <- trnL_fall16 %>% 
  full_join(trnL_spring17) %>% 
  full_join(trnL_fall17) %>% 
  full_join(vouchers)

# WHERE ARE MY PLANT VOUCHERS? #
# still missing ~ 27 plant vouchers
trnL_plant_vouchers <- trnL_props %>% 
  select(OTU, one_of(as.character(plant_vial_id$vial_id)))


