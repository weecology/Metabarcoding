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
# TODO: figure out what to do with the output - DONE
# TODO: Write script to only BLAST OTUs not already in
#       NoBlast_blastoutput.csv 
        
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
write.csv(clean.file, "NoBlast_blastoutput.csv")
