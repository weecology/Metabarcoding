# ITS2 plant samples OTUs
# EKB
# January 23, 2020

# LIBRARIES
library(tidyverse)

# DATA
trnL_reads <- read_csv("Data/SequencedData/Plants/ProcessedData/trnL_reads.csv")
trnL_totals <- read_csv("Data/SequencedData/Plants/ProcessedData/trnL_totals.csv")
sample_id <- read_csv("Data/CollectionData/vial_id.csv") %>% 
  filter(sample_type == "plant")
collection_id <- read_csv("Data/CollectionData/plant_voucher_collection.csv") %>% 
  select(vial_barcode, sci_name_profID) %>% 
  drop_na()

# Get only plant vouchers
trnL_reads <- filter(trnL_reads, SampleID %in% sample_id$vial_id)
trnL_totals <- filter(trnL_totals, SampleID %in% sample_id$vial_id)

# Prep for loop

voucher_ids <- unique(trnL_reads$SampleID)
#test <- voucher_ids[1:10]

OTU_list <- list() # can use plyr::ldply(list, data.frame) to merge into dataframe when done
p.list <- list()
list_num = 1
p.list_num = 1

for (i in 1:length(test)) {
  
  id <- voucher_ids[i]
  id_OTUs <- trnL_reads[(trnL_reads$SampleID == id), ]
  print(collection_id[collection_id$vial_barcode == id,])
  readline(prompt = "Press [enter] to continue")
  
  print(trnL_totals[trnL_totals$SampleID == id,])
  readline(prompt = "Press [enter] to continue")
  
  print(id_OTUs)
  hist(id_OTUs$Reads)
  
  print("Input the number of OTUs to keep [0-9] or type 'pass'")
  
  rows.to.keep = readline()
  if(rows.to.keep == 0) {
    # put NA in OTU_list
    OTU_list[[list_num]] <- data.frame(OTU = NA,
                                       SampleID = id,
                                       Reads = NA)
    list_num <- list_num + 1
  } else if (rows.to.keep %in% 1:20) {
    # put that number of rows in OTU_list
    OTU_list[[list_num]] <- id_OTUs[1:rows.to.keep,]
    list_num <- list_num + 1
  } else {
    # put that sample id in the "pass" list to revist
    p.list[[p.list_num]] <- id
    p.list_num <- p.list_num + 1
  }
  readline(prompt = "Press [enter] to continue")
}
