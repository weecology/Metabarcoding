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

# Fall 2017 
# good

#  Spring 2017 #
spring2017 <- spring2017[, -c(155:157)]

spring2017$Sum <- rowSums(spring2017[,-c(1:2)])

total_reads_row <- c("Total_Reads", NA, colSums(spring2017[,-c(1:2)]))
spring2017 <- rbind(spring2017, total_reads_row)

# Spring 2018 #
total_reads_row <- c("Total_Reads", NA, colSums(spring2018[,-c(1:2)]))
spring2018 <- rbind(spring2018, total_reads_row)

# Trap & Bait #
trapandbait$Sum <- rowSums(trapandbait[,-c(1:2)])

# PROCESS FILES #

df_list <- list(fall2017 = fall2017, 
                spring2017 = spring2017, 
                spring2018 = spring2018, 
                trapandbait = trapandbait)

for (i in names(df_list)){
  
  # fix naming (i.e., remove Wisely from vial IDs)
  colnames(df_list[[i]]) <- sub(".Wisely|.Bledsoe|Bledsoe.", "", colnames(df_list[[i]]))
  colnames(df_list[[i]]) <- sub("[.]2", "", colnames(df_list[[i]]))
  colnames(df_list[[i]]) <- sub("[.]1", "", colnames(df_list[[i]]))
  colnames(df_list[[i]])[1] <- "OTU"
  
  assign(i, df_list[[i]])
  
  write_csv(df_list[[i]], 
            paste0("Data/SequencedData/Plants/RawData/PreppedFiles/ITS2/", i, ".csv"))
  
}
