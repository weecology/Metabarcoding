# First Pass at Plotting 
# Fecal samples, ITS
# April 2020 EKB

# LIBRARIES, SOURCE CODE & DATA #

library(tidyverse)
library(vegan)
source("Scripts/functions.R")

reads <- read_csv("Data/SequencedData/Plants/ProcessedData/trnL_reads.csv")
totals <- read_csv("Data/SequencedData/Plants/ProcessedData/trnL_totals.csv")
samples <- read_csv("Data/CollectionData/fecal_sample_collection.csv")
