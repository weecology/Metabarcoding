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
# TODO: put OTUname and sequence together in a FASTA format - PARTIAL
#           current status: for millet (OTU1), I paste the info
#           together in a fasta consistent fashion that can be
#           fed to blastSequences(). Currently this is a one
#           sequence at a time affair. It needs to be automated
#           for multiple sequences
# TODO: figure out how to submit those sequences to GenBank - PARTIAL
#         SKME: currently it can submit the sequence for OTU1 and
#               provides the results in a data.frame. It IDs
#               Panicum millaceum, as we found using the BLAST GUI
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
#library(RCurl)
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
OTU1 = noblast_OTUs %>% filter(OTU_its == "OTU1")

### These packages seem to fight with dplyr, 
###   so I don't load them until I need them
source("https://bioconductor.org/biocLite.R")
biocLite("annotate")
library(annotate)

### Pastes everything together for a fasta format

header = paste(">",OTU1$OTU_its, sep="")
data = as.character(paste(header,OTU1$sequence_its,sep='\n'))

### Queries BLAST
output = blastSequences(x=data,
                        hitListSize = 20, as='data.frame')

