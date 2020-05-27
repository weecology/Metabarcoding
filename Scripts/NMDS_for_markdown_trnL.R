# NMDS trnL Markdown Plotting
# EKB
# May 2020

library(tidyverse)
library(vegan)
source('Scripts/functions.R')

reads <- read_csv("Data/SequencedData/Plants/ProcessedData/trnL_reads.csv")
totals <- read_csv("Data/SequencedData/Plants/ProcessedData/trnL_totals.csv")
samples <- read_csv("Data/CollectionData/fecal_sample_collection.csv")

# colorblind friendly palette
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# DATA PREP #

test_list_2017 <- filter_reads_data(samples, reads, totals, reads_min = 2000,
                                    yr = 2017, rel_reads_min = 0.05)
test_list2 <- data_prep_multivariate(test_list_2017)

# remove OTU150 (millet)
test_list2[[1]] <- test_list2[[1]] %>% select(-OTU150)

# remove outlier (S013067)
test_list2[[1]] <- test_list2$reads_spread[-105,]
test_list2[[2]] <- test_list2$sampleID[-105]
test_list2[[3]] <- test_list2$groups[-105,]

# make binary
test_list2[[1]] <- binarize(test_list2[[1]])

# RUN NMDS #
dist_trnL <- metaMDS(test_list2[[1]], distance = "bray", trymax = 250, k = 2)
# dist_trnL_converge <- run_metaMDS_til_converge(test_list2, dist_trnL, "bray", 2)
stressplot(dist_trnL)

# PLOT NMDS # 
NMDS_list <- NMDS_plotting_prep(test_list2, dist_trnL) 
plot_NMDS_ggplot2(NMDS_list)

