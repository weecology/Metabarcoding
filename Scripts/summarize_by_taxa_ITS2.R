# Summarize By Taxa - ITS2
# EKB
# May 2020

# LIBRARIES AND DATA # 

library(tidyverse)

taxa_data <- read_csv("Data/SequencedData/Plants/ProcessedData/ITS2_reads_with_lineage.csv")

# PULL OUT CONSENSUS LINEAGE #

taxa_data <- taxa_data %>% select(OTU, SampleID, Reads, DataFrame, ConsensusLineage)

# taxonomy dataframe

for(this_level in c('k','p','c','o','f','g','s')){
  
  # separate taxa into columns 
  step_one = sapply(strsplit(as.character(taxa_data$ConsensusLineage), paste0(this_level,'__')), 
                    '[', 2)
  step_two = sapply(strsplit(step_one, ';'), '[', 1)
  taxa_data[,this_level] = step_two
  
}

# rename columns
taxa_data_ITS2 <- rename(taxa_data, Domain = k, Clade1 = p, Class = c, 
                         Order = o, Family = f, Genus = g, Species = s) %>% 
  na_if("None")

# What do we do with the "unclassified" data? Remove it?