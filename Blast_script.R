#==========================================
# Blast_script.R
#
#   Code in this file extracts OTU sequences
#   from .fna file and runs them against
#   GenBank (or at least that's the goal)
#==========================================

################  TODO List: Description - STATUS, comments underneath

# TODO: read & extract OTUS from noblast file - DONE
# TODO: extract noblast OTU sequences from ITS .fna file - DONE
# TODO: put OTUname and sequence together in a FASTA format - DONE
# TODO: figure out how to submit those sequences to GenBank - DONE
# TODO: figure out what to do with the output - NOT EVEN CLOSE
#         SKME: I have no idea what any of that output means. I also
#               have no idea how we're going to extract the taxon name
#               because the species name is embedded with info on where
#               the species came from.
# TODO: extract 98% match taxonomies - MAYBE
#         SKME: I think blastSequence() can be used to specify match level
#               if we know what the match metric to use is.

#######################   Main Code

library(dplyr)
# to run this script, you will need to download bioconductor first:


### Reads no blast file, extracts the OTU.IDs
noblast = read.csv("./Plants/ITS_no_blast.csv", stringsAsFactors = FALSE)
OTUs = noblast %>% select(OTU.ID)

### Reads the .csv of the .fna file and extracts
###   only the sequences for the OTUs in the no
###   blast file

allITSseqs = read.csv("./Plants/ITS_from_fna.csv", stringsAsFactors = FALSE) 
noblast_OTUs = allITSseqs %>% filter(OTU_its %in% OTUs$OTU.ID)

### Before working with all 451, I'm testing out the code with
###  just OTU1, which we know is millet
OTU1 = noblast_OTUs %>% filter(OTU_its %in% c("OTU1","OTU101"))

### These packages seem to fight with dplyr, 
###   so I don't load them until I need them

source("https://bioconductor.org/biocLite.R")
biocLite("annotate")
library(annotate)

### Queries BLAST
###   Pastes OTU ID and sequence together for a fasta format
###   Submits to BLAST and records output
 
file = c()
num_seq = nrow(OTU1)
for(i in 1:num_seq){
  print(paste("Number of sequences remaining:",num_seq-(i-1),sep=" "))
  header = as.character(paste(">",OTU1$OTU_its[i], sep=""))
  data = as.character(paste(header,OTU1$sequence_its[i],sep="\n"))
  output = blastSequences(x=data,timeout = 220,
                          hitListSize = 20, as='data.frame')
  file = rbind(file,output)
  print(paste(OTU1$OTU_its[i], "complete", sep = " "))
  }
write.csv(file, "Blast_Results_ITS.csv")

