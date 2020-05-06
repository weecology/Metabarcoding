# Make LONG FORM ITS2 data
# EKB
# Jan 23, 2020

# LIBRARIES
library(tidyverse)
source("Scripts/functions.R")

# DATA
df_list <- read_in_ITS2_files()
names(df_list) <- c("fall2017", "spring2017", "spring2018", "trapandbait")

# Process Data

reads_list <- list()
total_list <- list()

for (i in 1:length(df_list)) {
  
  df <- df_list[[i]]
  
  df_total_reads <- as.vector(df[nrow(df),]) %>% 
    select(-ConsensusLineage) 
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
    pivot_longer(cols = c(3:ncol(df)), names_to = "SampleID", values_to = "Reads") %>% 
    filter(Reads > 0)
  df$DataFrame <- names(df_list[i])
  reads_list[[i]] <- df
  
}

ITS2_reads <- plyr::ldply(reads_list, data.frame) %>% 
  arrange(SampleID, (desc(Reads)))
ITS2_totals <- plyr::ldply(total_list, data.frame)

write_csv(ITS2_reads, "Data/SequencedData/Plants/ProcessedData/ITS2_reads_with_lineage.csv")
write_csv(ITS2_totals, "Data/SequencedData/Plants/ProcessedData/ITS2_totals.csv")
