# Update vial ID file

library(tidyverse)

# files
plants <- read_csv('Data/CollectionData/plant_voucher_collection.csv')
fecal <- read_csv('Data/CollectionData/fecal_sample_collection.csv')

# gather the data you need
plants <- plants %>% 
  select(vial_id = vial_barcode, sample_id = label_number)
plants$sample_type <- "plant"

fecal <- fecal %>% 
  select(vial_id = vial_barcode, sample_id = PIT_tag)
fecal$sample_type <- "fecal"

# merge together
vial_id_key <- rbind(plants, fecal) %>% 
  drop_na(vial_id)
write_csv(vial_id_key, "Data/CollectionData/vial_id.csv")
