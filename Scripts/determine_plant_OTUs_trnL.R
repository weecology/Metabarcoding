# trnL Plant Sample OTUs
# EKB
# January 22, 2020

# LIBRARIES
library(tidyverse)

# DATA
trnL_reads <- read_csv("Data/SequencedData/Plants/ProcessedData/trnL_reads.csv")
trnL_totals <- read_csv("Data/SequencedData/Plants/ProcessedData/trnL_totals.csv")
sample_id <- read_csv("Data/CollectionData/vial_id.csv") %>% 
  filter(sample_type == "plant")

# Get only plant vouchers
trnL_reads <- filter(trnL_reads, SampleID %in% sample_id$vial_id)
trnL_totals <- filter(trnL_totals, SampleID %in% sample_id$vial_id)

# Make Histograms

voucher_ids <- unique(trnL_reads$SampleID)

OTU_list <- list() # can you plyr::ldply(list, data.frame) to merge into dataframe when done
passed_list <- list()

for (i in 1:length(voucher_ids)) {
  
  id <- voucher_ids[i]
  id_OTUs <- trnL_reads[(trnL_reads$SampleID == id), ]
  print(id_OTUs)
  hist(id_OTUs$Reads)
  
  print("Input the number of OTUs to keep [0-9] or type 'pass'")
  
  rows.to.keep = readline()
  if(rows.to.keep == 0) {
    # put NA in dataframe
    # ## To remove a star:
    # ws[which(ws$tag == extrastar[i, 'tag']), 'note2'] <- NA
    # print(ws[which(ws$tag == extrastar[i, 'tag']), ])
    # print('Remember to record on datasheet + in notebook!')
  } else if (rows.to.keep %in% 1:20){
    # put those OTUs in the dataframe with that sample number (list)
  } else {
    # make list of 
  }
  readline(prompt = "Press [enter] to continue")
}

print("You've gone through all the vouchers. Type 'Y' if you want to revisit the ones you passed")
revisit.passed = readline()
if (revisit.passed == 'Y') {
  for (i in 1:length(passed_list))
}
