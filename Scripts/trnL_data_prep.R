# Prep trnL data from Jonah Inc.
# Ellen Bledsoe
# Jan 2017

########################
# LIBRARIES

library(dplyr)
library(stringr)

########################
# LOAD FILES

trnL <- read.csv("C:/Users/ellen.bledsoe/Dropbox (UFL)/Portal/PORTAL_primary_data/DNA/Results_Jonah/Plants/March2017/otu.table.trnL.cr99.exh.jv38.tax.csv", header = T, stringsAsFactors = F)
samples <- read.csv("CollectionData/fecal_sample_collection.csv", header = T, stringsAsFactors = F) 

########################
# CLEAN DATA

# taxonomy dataframe
taxa_trnL <- select(trnL, OTU_ID, ConsensusLineage) %>% 
  rename(OTU.ID = OTU_ID)

for(this_level in c('d','k','p','c','o','f','g','s')){
  # separate taxa into columns
  step_one=sapply(strsplit(as.character(taxa_trnL$ConsensusLineage), paste0(this_level,'__')), '[', 2)
  step_two=sapply(strsplit(step_one, ';'), '[', 1)
  taxa_trnL[,this_level]=step_two
}

taxa_trnL <- rename(taxa_trnL, Family = o, Genus = g, Species = s)

# reads dataframe

reads <- trnL[,c(1,8:54)] %>% rename(OTU.ID = OTU_ID)
reads <- tidyr::gather(reads, "Sample", "Reads", 2:48) %>% 
  filter(Reads != 0)

