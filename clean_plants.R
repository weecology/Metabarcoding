# Clean and Standardize Plant Names
# EKB
# July 2017

library(dplyr)

# read and prep master plant list
master <- read.csv("Portal_plant_master_list.csv", stringsAsFactors = F, strip.white = T)
master <- tidyr::unite(master, sci_name, Genus, Species, sep = " ", remove = F)

# read and prep plant collection file
collection <- read.csv("CollectionData/plant_voucher_collection.csv", stringsAsFactors = F, strip.white = T)
collection <- filter(collection, 
                     label_number != 'millet',
                     DNA == 'Y')
fillers <- collection$sci_name_fieldID[106:127]
collection$sci_name_profID[106:127] = fillers

### in master list, not collection

# by species code       - 58 species
in_master_spcode <- anti_join(master, collection, by = c("Sp.Code" = "sp_code"))

# by sciName            - 75 species
in_master_sciName <- anti_join(master, collection, by = c("sci_name" = "sci_name_profID"))

### in collection, not master

# by species code       - 9 species
in_coll_spcode <- anti_join(collection, master, by = c("sp_code" = "Sp.Code"))
# changed in collection:
  # aris tern -> aris hamu (but ternipes isn't accepted name?)
# unconfirmed: 
  # cymo mult
  # plan sp
  # chen sp
# TO BE added to master:
  # Setaria leucopila (seta leuc)
  # Cuscuta umbellata (cusc umbe)
  # Sida neomexicana (sida neom)
  # Abutilon parvulum (abut parv)
  # Delphinium wootonii (delp woot)

# by sciName            - 29 species
in_coll_sciName <- anti_join(collection, master, by = c("sci_name_profID" = "sci_name"))
# changed in collection:
  # Lappula redowskii -> Lappula occidentalis
  # Sphaeralcea coccinea -> Sphaeralcea hastulata
  # Salsola tragus -> Kali tragus
  # Talinum aura -> Phemeranthus aurantiacus
  # erio lemm -> erio acum
  # Bouteloua barbata -> Chondrosum barbatum
  # acac greg -> mimo acul
  # Dithyrea wislizeni -> Dimorphocarpa wislizeni
  # Dichelostemma pulchellum -> Dichelostemma capitatum
# TO BE added/changed in master:
  # Polygala tweedyi -> Polygala lindheimeri var. parvifolia
  # Erioneuron pulchellum -> Dasyochloa pulchella
  # Talinum aura & angu -> Phemeranthus aurantiacus
  # Tetraclea coulteri -> Clerodendrum coulteri
  # Eriochloa acuminata (erio acum)
  # Acacia constricta -> Vachellia constricta
  # mimo acul, Mimosa aculeaticarpa
# unconfirmed
  # Senegalia (Acacia) greggii vs. Mimosa aculeaticarpa