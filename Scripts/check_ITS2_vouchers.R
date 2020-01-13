# Check for ITS2 Data for All Samples
# Jan 13 2020
# EKB

# This script is to confirm that I have trnL data for all plant and fecal samples

# LIBRARIES #
library(tidyverse)

#==============================================================================#
# DATA PREP 

# Get Vial IDs #
vial_id <- read_csv("Data/CollectionData/vial_id.csv")
plant_vial_id <- vial_id %>% filter(sample_type == "plant")
fecal_vial_id <- vial_id %>% filter(sample_type == "fecal")

# Get ITS2 Data #

# Fall 2017 -- corrected OTUs 
fall2017 <- read_csv("Data/SequencedData/Plants/RawData/ITS2_Fall2017_correctedOTUs.csv")
fall2017 <- fall2017[,-c(2, length(fall2017))]
fall2017 <- fall2017[-c((nrow(fall2017)-1):nrow(fall2017)),]

# Spring 2017 -- corrected OTUs
spring2017 <- read_csv("Data/SequencedData/Plants/RawData/ITS2_Spring2017_correctedOTUs.csv")
spring2017 <- spring2017[,-(length(spring2017))]

# Spring 2018 
spring2018 <- read_csv("Data/SequencedData/Plants/RawData/ITS2_Spring2018.csv")
spring2018 <- spring2018[,-c(2, length(spring2018))]

# Trap & Bait Test
trap_bait <- read_csv("Data/SequencedData/Plants/RawData/ITS2_Trap_Bait_Test.csv")
trap_bait <- trap_bait[-(nrow(trap_bait)), -2]

# PROCESS FILES #

df_list <- list(fall2017 = fall2017, 
                spring2017 = spring2017, 
                spring2018 = spring2018,
                trap_bait = trap_bait)

for (i in names(df_list)){
  
  # fix naming (i.e., remove Wisely from vial IDs)
  colnames(df_list[[i]]) <- sub(".Wisely|.Bledsoe|Bledsoe.", "", colnames(df_list[[i]]))
  colnames(df_list[[i]]) <- sub("[.]1", "", colnames(df_list[[i]]))
  colnames(df_list[[i]])[1] <- "OTU"
  
  assign(i, df_list[[i]])
  
}

combined_df <- df_list[[1]] %>% 
  full_join(df_list[[2]]) %>% 
  full_join(df_list[[3]]) %>% 
  full_join(df_list[[4]])

#==============================================================================#
# CHECK FOR SAMPLES 

### Plant Vouchers ###
# all plant vouchers except millet and the missing spha hast are accounted for

plant_vouchers <- combined_df %>% 
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
setequal(meta_fecal$Sample_Barcode, fecal_vial_id$vial_id) 

# missing some samples
fecal_samples <- combined_df %>% 
  select(OTU, one_of(as.character(fecal_vial_id$vial_id)))
diff <- setdiff(meta_fecal$Sample_Barcode, colnames(fecal_samples[,-1]))

# looking for patterns in missing samples (I don't think all got run)
#   - seems like the ones from 460 that are missing are all the fresh ones
#     that we decided we didn't need based on the trap & bait tests (-460)
#   - missing sample from 454 and samples from 466 probably had moisture in them
#   - missing a couple additional (4, I think), but that's because some just didn't run
# GOOD TO GO
fecal_samples <- read_csv("Data/CollectionData/fecal_sample_collection.csv")
missing <- filter(fecal_samples, vial_barcode %in% diff)
