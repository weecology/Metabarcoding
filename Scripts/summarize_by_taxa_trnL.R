# Summarize By Taxa - trnL
# EKB
# May 2020

# LIBRARIES AND DATA # 

library(tidyverse)

taxa_data <- read_csv("Data/SequencedData/Plants/ProcessedData/trnL_reads_with_lineage.csv")

# PULL OUT CONSENSUS LINEAGE #

taxa_data <- taxa_data %>% select(OTU, SampleID, Reads, ConsensusLineage)

# taxonomy dataframe

for(this_level in c('d','k','p','c','o','f','g','s')){
  
  # separate taxa into columns 
  # this only works for some data (sep = '__')
  step_one = sapply(strsplit(as.character(taxa_data$ConsensusLineage), 
                             paste0(this_level,'__')), '[', 2)
  step_two = sapply(strsplit(step_one, ';'), '[', 1)
  taxa_data[,this_level] = step_two

}

# pull out the rows for which the for loop above didn't work

df1 <- taxa_data[!is.na(taxa_data$d),]
df2 <- taxa_data[is.na(taxa_data$d),]

for(this_level in c('d','k','p','c','o','f','g','s')){
  
  # separate taxa into columns 
  # this only works for some data (sep = ':')
  step_one = sapply(strsplit(as.character(df2$ConsensusLineage), 
                             paste0(this_level,':')), '[', 2)
  step_two = sapply(strsplit(step_one, ','), '[', 1)
  step_three = sapply(strsplit(step_two, ';'), '[', 1)
  df2[,this_level] = step_three

}

# row bind the two dataframes back together
taxa_data_trnL <- bind_rows(df1, df2)

# rename columns
taxa_data_trnL <- rename(taxa_data_trnL, Kingdom = d, Clade1 = k, Clade2 = p, 
                         Order = c, Family = o, Subfamily = f, Genus = g, 
                         Species = s) %>% 
  na_if("")

# Add WeeTUs for Summarizing
taxa_data_trnL <- taxa_data_trnL %>% 
  mutate(WTU.kingdom = group_indices(., Kingdom),
         WTU.clade1 = group_indices(., Kingdom, Clade1),
         WTU.clade2 = group_indices(., Kingdom, Clade1, Clade2),
         WTU.order = group_indices(., Kingdom, Clade1, Clade2, Order),
         WTU.family = group_indices(., Kingdom, Clade1, Clade2, Order, Family),
         WTU.subfamily = group_indices(., Kingdom, Clade1, Clade2, Order, 
                                       Family, Subfamily),
         WTU.genus = group_indices(., Kingdom, Clade1, Clade2, Order, Family, 
                                   Subfamily, Genus),
         WTU.species = group_indices(., Kingdom, Clade1, Clade2, Order, Family, 
                                     Subfamily, Genus, Species))

#write_csv(taxa_data_trnL, "Data/SequencedData/Plants/ProcessedData/trnL_reads_WeeTU.csv")
