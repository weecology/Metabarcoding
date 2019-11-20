# Prep trnL files for easy use
# 19 November 2019
# EKB

# This script is to make the trnL files that will go in the "PreppedFiles" 
# folder in Data/SequencedData/Plants/RawData. These CSVs will then be easy to 
# write functiosn to use on, as they will be a consistent structure

# LIBRARY #
library(tidyverse)

# READ FILES #
fall2016 <- read_csv("Data/SequencedData/Plants/RawData/trnL_Fall2016.csv")
fall2017 <- read_csv("Data/SequencedData/Plants/RawData/trnL_Fall2017.csv")
millet <- read_csv("Data/SequencedData/Plants/RawData/trnL_millet.csv")
refsamples <- read_csv("Data/SequencedData/Plants/RawData/trnL_refsamples.csv")
spring2017 <- read_csv("Data/SequencedData/Plants/RawData/trnL_Spring2017.csv")
spring2018 <- read_csv("Data/SequencedData/Plants/RawData/trnL_Spring2018.csv")
trapandbait <- read_csv("Data/SequencedData/Plants/RawData/trnL_trapandbait.csv")

# PREP FILES #

# Fall 2016 #
# convert fall2016 to reads to match others
fall2016 <- fall2016[,-c(2:6)]
fall2016[nrow(fall2016), ncol(fall2016)] <- rowSums(fall2016[nrow(fall2016), -c(1:2, ncol(fall2016))])
fall2016_totalreads <- fall2016[nrow(fall2016),] 
fall2016_OTU <- fall2016[-nrow(fall2016), c(1:2)]

fall2016 <- mapply('*', fall2016[-nrow(fall2016), -c(1:2)], fall2016_totalreads[-c(1:2)]) %>% 
  round(digits = 0)
fall2016 <- cbind(fall2016_OTU, fall2016) %>% 
  rbind(fall2016_totalreads) %>% 
  as_tibble()
fall2016[nrow(fall2016), ncol(fall2016)] <- NA

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

# PROCESS FILES #

df_list <- list(fall2016 = fall2016, 
                fall2017 = fall2017, 
                millet = millet, 
                refsamples = refsamples, 
                spring2017 = spring2017, 
                spring2018 = spring2018, 
                trapandbait = trapandbait)

for (i in names(df_list)){
  
  # fix naming (i.e., remove Wisely from vial IDs)
  colnames(df_list[[i]]) <- sub(".Wisely|.Bledsoe|Bledsoe.", "", colnames(df_list[[i]]))
  colnames(df_list[[i]]) <- sub("[.]2", "", colnames(df_list[[i]]))
  colnames(df_list[[i]])[1] <- "OTU"
  
  assign(i, df_list[[i]])
  
  write_csv(df_list[[i]], 
            paste0("Data/SequencedData/Plants/RawData/PreppedFiles/trnL/", i, ".csv"))
  
}
