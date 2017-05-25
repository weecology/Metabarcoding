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

joined <- select(joined, -ConsensusLineage) %>%
  rename(Family = f, Genus = g, SciName = s) %>% 
  filter(k != "Viruses")

#########################################
# MATCH OTU.ID w/ SEQUENCE

seq <- read.csv("./SequencedData/Plants/ITS_sequences_from_fna.csv", header = T, stringsAsFactors = F)
seq <- rename(seq, OTU.ID = OTU_its)
OTU_for_dups <- joined$OTU.ID
seq_to_BLAST <- filter(seq, OTU.ID %in% OTU_for_dups)

#########################################
# RUN THROUGH BLAST

completed_blasts = read.csv("./SequencedData/Plants/reBlast_blastoutput.csv", stringsAsFactors = FALSE)
completed_OTUs = unique(completed_blasts$OTU.ID)
OTUs_forBLAST = noblast_OTUs %>% filter(!(OTU_its %in% completed_OTUs))

### These packages seem to fight with dplyr, 
###   so I don't load them until I need them


source("https://bioconductor.org/biocLite.R")
biocLite("annotate")
library(annotate)

### Queries BLAST
###   Pastes OTU ID and sequence together for a fasta format
###   Submits to BLAST and records output

file = c()
num_seq = nrow(OTUs_forBLAST)
for(i in 1:num_seq){
  print(paste("Number of sequences remaining:",num_seq-(i-1),sep=" "))
  header = as.character(paste(">",OTUs_forBLAST$OTU_its[i], sep=""))
  data = as.character(paste(header,OTUs_forBLAST$sequence_its[i],sep="\n"))
  output = blastSequences(x=data,timeout = 220,
                          hitListSize = 20, as='data.frame')
  file = rbind(file,output)
  print(paste(OTUs_forBLAST$OTU_its[i], "complete", sep = " "))
}

### Formats and Write table from BLAST
### Selects & formats only relevant columns
### Calculates Identity% & Query Coverage

names(file)[3] = 'OTU.ID'
names(file)[4] = 'Query.length'
names(file)[25] = 'Hsp.length'

clean.file = dplyr::select(file, OTU.ID, Query.length, Hit_id, Hit_def,
                           Hit_len, Hsp_evalue, 
                           Hsp_identity, Hsp.length)
clean.file = clean.file %>% mutate_each_(funs(as.integer), c("Query.length","Hit_len",
                                                             "Hsp_identity","Hsp.length"))
clean.file = clean.file %>% mutate_each_(funs(as.numeric), "Hsp_evalue")

clean.file = clean.file %>% mutate(identity_percent = 100* (Hsp_identity/Hsp.length),
                                   Query.cover = 100* (Hsp.length/Query.length))

completed_blasts = rbind(completed_blasts,clean.file)

write.csv(completed_blasts, "./Plants/NoBlast_blastoutput.csv", row.names=FALSE)
