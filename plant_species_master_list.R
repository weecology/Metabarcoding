# Master Plant List and trnL/ITS Stats
# EKB
# May 2017

# LIBRARIES
library(dplyr)

# FILES

plant_list <- read.csv("./Portal_plant_master_list.csv")
plant_collection <- read.csv("./CollectionData/plant_voucher_collection.csv")

# PORTAL PLANTS MASTER LIST

unknowns <- filter(plant_list, Species %in% c('sp.', 'spp.', 'sp. 2') | Genus %in% c('unknown', 'Unknown'))

plant_list <- anti_join(plant_list, unknowns)
