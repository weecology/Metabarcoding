# Check for trnL Data for All Samples
# 17 November 2019
# EKB

# This script is to confirm that I have trnL data for all plant and fecal samples

# LIBRARIES #
library(tidyverse)

#==============================================================================#
# DATA PREP 

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
# from the Portal Dropbox "reference_samples" folder

vouchers <- read_csv("Data/SequencedData/Plants/RawData/trnL_refsamples.csv")

# remove extra columns
vouchers <- vouchers[, -2]
colnames(vouchers)[1] <- "OTU"
colnames(vouchers) <- sub(".Wisely", "", colnames(vouchers))

vouchers$Sum <- rowSums(vouchers[,-1])
vouchers_OTUsum <- vouchers[, c(1, length(vouchers))]
vouchers_total_reads <- vouchers[nrow(vouchers),]
vouchers <- vouchers[-nrow(vouchers),-length(vouchers)]

## trnL Trap & Bait Test ##
# from the Portal Dropbox Trap_Bait_Test folder

trnL_trapbait <- read_csv("Data/SequencedData/Plants/RawData/trnL_trapandbait.csv")

# fix naming (i.e., remove Wisely from vial IDs)
colnames(trnL_trapbait) <- sub(".Wisely|.Bledsoe", "", colnames(trnL_trapbait))
colnames(trnL_trapbait)[1] <- "OTU"

# remove extra columns and rows
trnL_trapbait <- trnL_trapbait[,-2]
trnL_trapbait$Sum <- rowSums(trnL_trapbait[,-1])
trnL_trapbait_OTUsum <- trnL_trapbait[, c(1, length(trnL_trapbait))]
trnL_trapbait_total_reads <- trnL_trapbait[nrow(trnL_trapbait),]
trnL_trapbait <- trnL_trapbait[-nrow(trnL_trapbait),-length(trnL_trapbait)]

## trnL Millet ##
# also from the Portal Dropbox Trap_Bait_Test folder

trnL_millet <- read_csv("Data/SequencedData/Plants/RawData/trnL_millet.csv")

# fix naming (i.e., remove Wisely from vial IDs)
colnames(trnL_millet) <- sub(".Wisely|.Bledsoe", "", colnames(trnL_millet))
colnames(trnL_millet)[1] <- "OTU"

# remove extra columns and rows
trnL_millet <- trnL_millet[,-2]
trnL_millet_total_reads <- trnL_millet[nrow(trnL_millet),]
trnL_millet <- trnL_millet[-nrow(trnL_millet),]

## trnL Spring 2018 ##

trnL_spring18 <- read_csv("Data/SequencedData/Plants/RawData/trnL_Spring2018.csv")

# fix naming (i.e., remove Wisely from vial IDs)
colnames(trnL_spring18) <- sub(".Wisely|.Bledsoe|Bledsoe.", "", colnames(trnL_spring18))
colnames(trnL_spring18) <- sub("[.]2", "", colnames(trnL_spring18))
colnames(trnL_spring18)[1] <- "OTU"

# remove extra columns and rows
trnL_spring18 <- trnL_spring18[,-2]
trnL_spring18 <- rbind(trnL_spring18, c("Total_Reads", colSums(trnL_spring18[,-1])))
trnL_spring18 <- trnL_spring18[!(trnL_spring18$Sum == 0),]
trnL_spring18_total_reads <- trnL_spring18[nrow(trnL_spring18),]
trnL_spring18_OTUsum <- trnL_spring18[, c(1, length(trnL_spring18))]
trnL_spring18 <- trnL_spring18[-(nrow(trnL_spring18)), -length(trnL_spring18)]


#==============================================================================#
# CHECK FOR SAMPLES 

# Note: these datasets are a mixture of reads and proportions
# I am combining them solely for the purpose of checking if I have all of them

# combine dataframes
trnL_props <- trnL_fall16 %>% 
  full_join(trnL_spring17) %>% 
  full_join(trnL_fall17) %>% 
  full_join(vouchers) %>% 
  full_join(trnL_trapbait) %>% 
  full_join(trnL_millet) %>% 
  full_join(trnL_spring18)

### Plant Vouchers ###
# all plant vouchers (except one? have been found) 
#  - S009904 (spha hast) is missing, but we have another spha hast so we're good
trnL_plant_vouchers <- trnL_props %>% 
  select(OTU, one_of(as.character(plant_vial_id$vial_id)))

### Fecal Samples ###

# not all fecal samples were sent in, so need to read in metadata sent to JV
metadata1 <- read_csv("C:/Users/ellen.bledsoe/Dropbox (UFL)/Portal/PORTAL_primary_data/DNA/metadata_jonah_2016.csv")
metadata2 <- read_csv("C:/Users/ellen.bledsoe/Dropbox (UFL)/Portal/PORTAL_primary_data/DNA/metadata_jonah_20170425.csv")
metadata3 <- read_csv("C:/Users/ellen.bledsoe/Dropbox (UFL)/Portal/PORTAL_primary_data/DNA/metadata_jonah_20171029.csv")
metadata4 <- read_csv("C:/Users/ellen.bledsoe/Dropbox (UFL)/Portal/PORTAL_primary_data/DNA/metadata_jonah_20180412.csv")

meta_meta <- metadata1 %>% 
  bind_rows(metadata2) %>% 
  bind_rows(metadata3) %>% 
  bind_rows(metadata4) %>% 
  rename("Sample_Type" = "Sample Type")

meta_fecal <- meta_meta %>% 
  filter(Sample_Type == "fecal")
setequal(meta_fecal$Sample_Barcode, fecal_vial_id$vial_id) # TRUE

# missing some samples
trnL_fecal_samples <- trnL_props %>% 
  select(OTU, one_of(as.character(fecal_vial_id$vial_id)))
trnL_diff <- setdiff(meta_fecal$Sample_Barcode, colnames(trnL_fecal_samples[,-1]))

# looking for patterns in missing samples (I don't think all got run)
#   - seems like the ones from 460 that are missing are all the fresh ones
#     that we decided we didn't need based on the trap & bait tests (-460)
#   - missing sample from 454 and samples from 466 probably had moisture in them
# SHOULD BE GOOD TO GO HERE
fecal_samples <- read_csv("Data/CollectionData/fecal_sample_collection.csv")
trnL_missing <- filter(fecal_samples, vial_barcode %in% trnL_diff)
