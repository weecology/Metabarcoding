# Make LONG FORM ITS2 data
# EKB
# Jan 23, 2020

# LIBRARIES
library(tidyverse)
source("Scripts/functions.R")

# DATA
df_list <- read_in_trnL_files()

# Process Data

reads_list <- list()
total_list <- list()

for (i in 1:length(df_list)) {
  
  df <- df_list[[i]]
  df <- df %>% select(-ConsensusLineage) 
  
  df_total_reads <- as.vector(df[nrow(df),])
  df2 <- data.frame(t(df_total_reads[-1]))
  colnames(df2) <- df_total_reads[, 1]
  df2 <- rownames_to_column(df2, var = "SampleID") %>% 
    filter(SampleID != "Sum")
  total_list[[i]] <- df2
  
  df <- df[-nrow(df),]
  df <- df %>% 
    filter(Sum > 0) %>% 
    select(-Sum)
  df <-  df %>% 
    pivot_longer(cols = c(2:ncol(df)), names_to = "SampleID", values_to = "Reads") %>% 
    filter(Reads > 0)
  reads_list[[i]] <- df
  
}

trnL_reads <- plyr::ldply(reads_list, data.frame) %>% 
  arrange(SampleID, (desc(Reads)))
trnL_totals <- plyr::ldply(total_list, data.frame)

write_csv(trnL_reads, "Data/SequencedData/Plants/ProcessedData/trnL_reads.csv")
write_csv(trnL_totals, "Data/SequencedData/Plants/ProcessedData/trnL_totals.csv")
