library(dplyr)
library(tidyr)


data = read.csv("TRNL_portalplants_consensusID.csv")
data = data[-c((nrow(data)-1),nrow(data)),]
data = data %>% 
  select(-Sum) %>% 
  gather(SampleID, TrnL.Reads, S008809.Wisely:S009099.Wisely)
data = data %>% separate(SampleID, c("Sample", "Name"))
data = data %>% separate(ConsensusLineage, c("domain", 
                                             "kingdom",
                                             "phylum", 
                                             "class", 
                                             "order",
                                             "family",
                                             "genus",
                                             "species"),
                         sep=";")
data$domain = substring(data$domain,4)
data$kingdom = substring(data$kingdom,5)
data$phylum = substring(data$phylum,5)
data$class = substring(data$class,5)
data$order = substring(data$order,5)
data$family = substring(data$family,5)
data$genus = substring(data$genus,5)
data$species = substring(data$species,5)

best_match = data %>% group_by(Sample) %>% filter(TrnL.Reads == max(TrnL.Reads))
