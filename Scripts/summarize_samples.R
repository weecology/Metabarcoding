library(tidyverse)

plant_vouchers <- read_csv("Data/CollectionData/plant_voucher_collection.csv")
fecal_samples <- read_csv("Data/CollectionData/fecal_sample_collection.csv")
portal_plants <- read_csv("Data/CollectionData/Portal_plant_species.csv")

colnames(plant_vouchers)
n_distinct(plant_vouchers$sci_name_profID)
n_distinct(plant_vouchers$sci_name_profID[plant_vouchers$DNA == 'Y'])

sp.code_portal <- plant_list$speciescode
sp.code_collection <- plant_vouchers$sp_code

setdiff(sp.code_collection, sp.code_portal)
setdiff(sp.code_portal, sp.code_collection)

fecal_samples <- fecal_samples %>% 
  filter(is.na(notes), sample_type == 'trap') %>% 
  group_by(species, plot) %>% 
  count()

sum(fecal_samples$n[fecal_samples$species %in% c("DM", "DO")])
sum(fecal_samples$n[fecal_samples$species == 'PP'])

fecal_samples %>% 
  filter(species == 'PP', plot %out% c(3, 15, 19, 20, 21, 22))
