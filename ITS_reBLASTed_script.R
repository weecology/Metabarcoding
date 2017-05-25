# ITS duplicates through BLAST
# Ellen Bledsoe
# May 2017

#########################################
# LIBRARIES

library(dplyr)

#########################################
# LOAD FILES

blast <- read.csv("./Plants/ITS_blast.csv", header = TRUE, na.strings = "")

#########################################
# GET UNIQUE OTUs for DUPLICATES

# Blast and no blast files

blast <- select(blast, -Sum)
blast <- blast[-c(463:464),]

# Pull out ConsensusLineage

taxa_its <- select(blast, OTU.ID, ConsensusLineage)

# make dataframe with only duplicate lineages
dups <- taxa_its$ConsensusLineage
dups <- as.data.frame(dups[duplicated(dups)])
colnames(dups) <- "ConsensusLineage"

# join with OTU.ID to get all unique combinations
joined <- semi_join(taxa_its, dups)

for(this_level in c('k','p','c','o','f','g','s')){
  # separate taxa into columns
  step_one=sapply(strsplit(as.character(joined$ConsensusLineage), paste0(this_level,'__')), '[', 2)
  step_two=sapply(strsplit(step_one, ';'), '[', 1)
  joined[,this_level]=step_two
}

joined <- select(taxa_its, -ConsensusLineage) %>%
  rename(Family = f, Genus = g, SciName = s)
