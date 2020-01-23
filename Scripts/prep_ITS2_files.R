# Prep ITS2 files for easy use
# 23 January 2020
# EKB

# This script is to make the ITS2 files that will go in the "PreppedFiles" 
# folder in Data/SequencedData/Plants/RawData. These CSVs will then be easy to 
# write functiosn to use on, as they will be a consistent structure

# LIBRARY #
library(tidyverse)

# READ FILES #
fall2017 <- read_csv("Data/SequencedData/Plants/RawData/ITS2_Fall2017_correctedOTUs.csv")
spring2017 <- read_csv("Data/SequencedData/Plants/RawData/ITS2_Spring2017_correctedOTUs.csv")
spring2018 <- read_csv("Data/SequencedData/Plants/RawData/ITS2_Spring2018.csv")
trapandbait <- read_csv("Data/SequencedData/Plants/RawData/ITS2_Trap_Bait_Test.csv")

# PREP FILES #

# Fall 2017 #
fall2017 <- fall2017[-(nrow(fall2017)-1),]
fall2017[nrow(fall2017),1] <- "Total_Reads"

# Millet #
millet$Sum <- millet$S009063.Wisely

# Reference Samples #
refsamples$Sum <- rowSums(refsamples[,-c(1:2)])

#  Spring 2017 #
spring2017$Sum <- rowSums(spring2017[,-c(1:2)])

# Spring 2018 #
total_reads_row <- c("Total_Reads", NA, colSums(spring2018[,-c(1:2)]))
spring2018 <- rbind(spring2018, total_reads_row)

# Trap & Bait #
trapandbait$Sum <- rowSums(trapandbait[,-c(1:2)])