# Make trnL SEED only file 
# 12 November 2019
# EKB

# PREP MATERIALS #

# Prep trnL reference csv into functional dataframe #

trnL <- readr::read_csv("Data/SequencedData/Plants/RawData/trnL_ClosedRef_JV13.csv")
colnames(trnL) <- c("OTU", "Seed", "ConsensusLineage")

# split ConsesusLineage column into usable columns
for(this_level in c('d','k','p','c','o','f','g','s')){
  # separate taxa into columns
  step_one <- sapply(strsplit(as.character(trnL$ConsensusLineage), paste0(this_level,':')), '[', 2)
  step_two <- sub(";", "", step_one)
  step_three <- sapply(strsplit(step_two, ','), '[', 1)
  trnL[,this_level] <- step_three
}

# split species column into two columns and filter for SEED rows
trnL_seeds <- trnL %>%
  dplyr::select(., -g) %>% 
  tidyr::separate(s, c("g", "s"), " ") %>% 
  dplyr::filter(Seed == "SEED")

# write to CSV
# readr::write_csv(trnL_seeds, "Data/SequencedData/Plants/RawData/trnL_SEED_only.csv")