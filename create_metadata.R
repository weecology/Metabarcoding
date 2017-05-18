# Get Metadata for Jonah Ventures
# EKB
# 4/25/2016

#============================================================

# LOAD LIBRARIES and DATA

library(dplyr)

meta <- read.csv("C:/Users/ellen.bledsoe/Dropbox/Portal/PORTAL_primary_data/DNA/metadata_jonah_20170425.csv", header = TRUE)
plant <- read.csv("C:/Users/ellen.bledsoe/Dropbox/Portal/PORTAL_primary_data/DNA/Plants/plant_collection.csv", header = TRUE, na.strings = "")
fecal <- read.csv("C:/Users/ellen.bledsoe/Dropbox/Portal/PORTAL_primary_data/DNA/Rodents/fecal_samples.csv", header = TRUE)

#============================================================

# PLANT DATA

plant <- plant %>%
  filter(year == 2017, label_number != 'millet') %>%
  tidyr::unite(Date.Collected, month, day, year, sep = "/") %>% 
  tidyr::drop_na(vial_barcode)

plant <- rename(plant, 
                Species = sci_name_fieldID,
                Sample.ID.Barcode = vial_barcode)
plant <- select(plant, Date.Collected, Sample.ID.Barcode,Species)

meta$Date.Collected = as.character(meta$Date.Collected)
meta$Sample.ID.Barcode = as.factor(meta$Sample.ID.Barcode)
meta$Species = as.factor(meta$Species)

meta <- left_join(plant, meta)

meta$Sample.Type = 'plant'

# FECAL DATA

fecal <- fecal %>% 
  filter(year == 2017) %>% 
  tidyr::unite(Date.Collected, month, day, year, sep = "/") %>% 
  rename(Sample.ID.Barcode = vial_barcode)

fecal$Species = NA
for (i in 1:length(fecal$species)){
  if (fecal$species[i] == 'DO'){
    fecal$Species[i] = 'Dipodomys ordii'
  } else if (fecal$species[i] == 'DM'){
    fecal$Species[i] = 'Dipodomys merriami'
  } else {
    fecal$Species[i] = 'Chaetodipus pencillatus'
  }
}

fecal <- select(fecal, Date.Collected, Sample.ID.Barcode, Species)

meta <- bind_rows(meta, fecal)
meta$Sample.Type[23:166] = 'fecal'

# OTHER COLUMNS

meta$Sample.Num <- 1:nrow(meta)
meta$Latitude = '31.93785'
meta$Longitude = '-109.08305'
meta$City = 'Portal'
meta$State = 'AZ'
meta$Country = 'USA'
meta$Habitat = 'High desert'

meta <- select(meta, Sample.Num, Sample.ID.Barcode, OwnerID:Sample.Type, Date.Collected, Latitude:Country, Species, Habitat:UserRef)
write.csv(meta, "C:/Users/ellen.bledsoe/Desktop/metadata_jonah_2017.csv")
