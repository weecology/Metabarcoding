# Plot data histograms to determine cutoffs
# 26 Feb 2020
# EKB

# LIBRARIES #
library(tidyverse)

#==============================================================================#
# DATA PREP #

# Vial IDs 
vial_id <- read_csv("Data/CollectionData/vial_id.csv")
plant_vial_id <- vial_id %>% filter(sample_type == "plant")
fecal_vial_id <- vial_id %>% filter(sample_type == "fecal")

# sequenced data
trnL_data <- read_csv("Data/SequencedData/Plants/ProcessedData/trnL_reads.csv")
ITS_data <- read_csv("Data/SequencedData/Plants/ProcessedData/ITS2_reads.csv")

trnL_totals <- read_csv("Data/SequencedData/Plants/ProcessedData/trnL_totals.csv")
ITS2_totals <- read_csv("Data/SequencedData/Plants/ProcessedData/ITS2_totals.csv")

trnL_data <- full_join(trnL_data, trnL_totals) %>% 
  mutate(Rel_Reads = (Reads/Total_Reads)*100)

#==============================================================================#
# SELECT AND PLOT RANDOM DATA #

# get 16 random samples
unique_fecal <- unique(fecal_vial_id$vial_id)
random16 <- sample(unique_fecal, 16) # 16 random samples for plotting

# filter trnL and ITS for only random samples
trnL_16 <- trnL_data %>% filter(SampleID %in% random16)
ITS_16 <- ITS_data %>% filter(SampleID %in% random16)


# plot histograms
ggplot(data = (trnL_data %>% filter(SampleID %in% sample(unique_fecal, 1)))) +
  geom_histogram(aes(x = Rel_Reads)) +
  scale_y_log10() +
  scale_x_log10()

ggplot(trnL_16) +
  geom_histogram(aes(x = Rel_Reads), binwidth = .1) +
  facet_wrap(~ SampleID, nrow = 4, scales = "free") +
  scale_y_log10() +
  scale_x_log10()
