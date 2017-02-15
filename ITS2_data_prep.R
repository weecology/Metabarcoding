# Prep ITS2 data from Jonah Inc.
# Ellen Bledsoe
# Jan 2017

########################
# LIBRARIES

library(dplyr)
library(ggplot2)

########################
# LOAD FILES

blast <- read.csv("C:/Users/ellen.bledsoe/Desktop/Git/Metagenomics/Plants/ITS_blast.csv", header = TRUE)
no_blast <- read.csv("C:/Users/ellen.bledsoe/Desktop/Git/Metagenomics/Plants/ITS_no_blast.csv", header = TRUE)
plants <- read.csv("C:/Users/ellen.bledsoe/Desktop/Git/Metagenomics/CollectionData/plant_voucher_collection.csv", header = TRUE)

########################
# CLEAN DATA

# Blast
blast <- select(blast, -Sum)
reads_b <- tidyr::gather(blast, "Sample", "Reads", S008809.Wisely:S009099.Wisely) %>%
           filter(Reads != 0)
reads_b <- tidyr::separate(reads_b, Sample, into = c("Sample", "Wisely")) %>%
           select(-Wisely)

# No Blast
no_blast <- select(no_blast, -Sum)
reads_nb <- tidyr::gather(no_blast, "Sample", "Reads", S008809.Wisely:S009099.Wisely) %>%
  filter(Reads != 0)
reads_nb <- tidyr::separate(reads_nb, Sample, into = c("Sample", "Wisely")) %>%
  select(-Wisely)

# vouchers dataframe

vouchers <- select(plants, vial_barcode, label_number) %>% 
            rename(Sample = vial_barcode)
vouchers_b <- semi_join(reads_b, vouchers, by = "Sample")
vouchers_nb <- semi_join(reads_nb, vouchers, by = "Sample")

# fecal sample dataframe

fecal_b <- anti_join(reads_b, vouchers, by = "Sample")
fecal_nb <- anti_join(reads_nb, vouchers, by = "Sample")

########################
# UNIQUE OTUs for VOUCHERS

blast_match <- vouchers_b %>% group_by(Sample) %>% filter(Reads == max(Reads))
no_blast_match <- vouchers_nb %>% group_by(Sample) %>% filter(Reads == max(Reads))

for(this_level in c('k','p','c','o','f','g','s')){
  # separate taxa into columns
  step_one=sapply(strsplit(as.character(blast_match$ConsensusLineage), paste0(this_level,'__')), '[', 2)
  step_two=sapply(strsplit(step_one, ';'), '[', 1)
  blast_match[,this_level]=step_two
}

length(unique(blast_match$s))
###########################
# WORK AREA


