# First Pass at Plotting 
# Fecal samples, trnL
# March 2020 EKB

# LIBRARIES, SOURCE CODE & DATA #

library(tidyverse)
library(vegan)
source("Scripts/functions.R")

reads <- read_csv("Data/SequencedData/Plants/ProcessedData/trnL_reads.csv")
totals <- read_csv("Data/SequencedData/Plants/ProcessedData/trnL_totals.csv")
samples <- read_csv("Data/CollectionData/fecal_sample_collection.csv")

# DATA PREP #

# add plot type to fecal collection data
# add group for plotting
# and remove samples that were part of the trap/bait test
samples <- add_plot_type(samples) %>% 
  add_plotting_group() %>% 
  filter(is.na(notes))

# select only fecal samples
fecal_id <- samples$vial_barcode

# add totals to reads df, filter out fecal samples and small totals
reads <- full_join(reads, totals)
reads <- reads %>% 
  filter(SampleID %in% fecal_id) %>% 
  filter(Total_Reads > 2000) %>% 
  mutate(Rel_Reads = Reads/Total_Reads)



