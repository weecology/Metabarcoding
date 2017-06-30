# Master Plant List and trnL/ITS Stats
# EKB
# May 2017

# LIBRARIES
library(RCurl)
library(dplyr)

# FILES

plant_list <- read.csv(text = getURL("https://raw.githubusercontent.com/weecology/PortalData/master/Plants/Portal_plant_species.csv"))
plant_collection <- read.csv(text = getURL("https://raw.githubusercontent.com/weecology/Metagenomics/master/CollectionData/plant_voucher_collection.csv"))

# PORTAL PLANTS MASTER LIST

# remove any unknowns
unknowns <- filter(plant_list, Species %in% c('sp.', 'spp.', 'sp. 2') | Genus %in% c('unknown', 'Unknown'))
plant_list <- anti_join(plant_list, unknowns) %>% 
  arrange(Sp.Code)

# remove any species listed twice
duplicates <- select(plant_list, Sp.Code) %>%
  group_by(Sp.Code) %>%
  summarise(count = n()) %>%
  filter(count > 1)
duplicates <- filter(plant_list, Sp.Code %in% c('xant spin', 'yucc elat', 'zinn gran', 'zinn pumi') & Community == 'Perennials')
plant_list <- anti_join(plant_list, duplicates)

for (i in 1:length(plant_list$Sp.Code)) {
  if (plant_list$Sp.Code[i] == 'zinn pumi'){
    plant_list$Species[i] == 'acerosa' & plant_list$Alt.Species[i] == 'pumila'
  } 
}  

#write.csv(plant_list, "./Portal_plant_master_list.csv", row.names = F)
