# ITS2 plant samples OTUs
# EKB
# January 23, 2020

# LIBRARIES
library(tidyverse)

# DATA
ITS2_reads <- read_csv("Data/SequencedData/Plants/ProcessedData/ITS2_reads.csv")
ITS2_totals <- read_csv("Data/SequencedData/Plants/ProcessedData/ITS2_totals.csv")
sample_id <- read_csv("Data/CollectionData/vial_id.csv") %>% 
  filter(sample_type == "plant")
collection_id <- read_csv("Data/CollectionData/plant_voucher_collection.csv") %>% 
  select(vial_barcode, sci_name_profID) %>% 
  drop_na()

# Get only plant vouchers
ITS2_reads <- filter(ITS2_reads, SampleID %in% sample_id$vial_id)
ITS2_totals <- filter(ITS2_totals, SampleID %in% sample_id$vial_id)

# Prep for loop

voucher_ids <- unique(ITS2_reads$SampleID)
test <- voucher_ids[1:10]

OTU_list <- list() # can use plyr::ldply(list, data.frame) to merge into dataframe when done
p.list <- list()
list_num = 1
p.list_num = 1

for (i in 1:length(voucher_ids)) {
  
  id <- voucher_ids[i]
  id_OTUs <- ITS2_reads[(ITS2_reads$SampleID == id), ]
  print(collection_id[collection_id$vial_barcode == id,])
  readline(prompt = "Press [enter] to continue")
  
  print(ITS2_totals[ITS2_totals$SampleID == id,])
  readline(prompt = "Press [enter] to continue")
  
  print(id_OTUs)
  hist(id_OTUs$Reads)
  
  print("Input the number of OTUs to keep [0-9] or type 'pass'")
  
  rows.to.keep = readline()
  if(rows.to.keep == 0) {
    # put NA in OTU_list
    OTU_list[[list_num]] <- data.frame(OTU = NA,
                                       SampleID = id,
                                       Reads = NA,
                                       DataFrame = NA)
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

# If you DO NOT want to revisit the passed IDs, run lines 69-73.
# If you DO want to go back to through the passed IDs, go to line 75

OTU_df <- plyr::ldply(OTU_list, data.frame)
passed_df <- plyr::ldply(p.list, data.frame)

write_csv(OTU_df, "Data/SequencedData/Plants/ProcessedData/ITS2_OTUs.csv")
write_csv(passed_df, "Data/SequencedData/Plants/ProcessedData/ITS2_passed.csv")

# To revisit the passed IDs

p.list.2 <- list()
p.list_num = 1
  
for (i in 1:length(p.list)) {
  
  id <- p.list[i]
  id_OTUs <- ITS2_reads[(ITS2_reads$SampleID == id), ]
  print(collection_id[collection_id$vial_barcode == id,])
  readline(prompt = "Press [enter] to continue")
  
  print(ITS2_totals[ITS2_totals$SampleID == id,])
  readline(prompt = "Press [enter] to continue")
  
  print(id_OTUs)
  hist(id_OTUs$Reads)
  
  print("Input the number of OTUs to keep [0-9] or type 'pass'")
  
  rows.to.keep = readline()
  if(rows.to.keep == 0) {
    # put NA in OTU_list
    OTU_list[[list_num]] <- data.frame(OTU = NA,
                                       SampleID = id,
                                       Reads = NA,
                                       DataFrame = NA)
    list_num <- list_num + 1
  } else if (rows.to.keep %in% 1:20) {
    # put that number of rows in OTU_list
    OTU_list[[list_num]] <- id_OTUs[1:rows.to.keep,]
    list_num <- list_num + 1
  } else {
    # put that sample id in the "pass" list to revist
    p.list.2[[p.list_num]] <- id
    p.list_num <- p.list_num + 1
  }
  readline(prompt = "Press [enter] to continue")
}

OTU_df <- plyr::ldply(OTU_list, data.frame)
passed_df <- plyr::ldply(p.list, data.frame)

write_csv(OTU_df, "Data/SequencedData/Plants/ProcessedData/ITS2_OTUs.csv")
write_csv(passed_df, "Data/SequencedData/Plants/ProcessedData/ITS2_passed.csv")






