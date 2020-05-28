# NMDS trnL Markdown Plotting
# EKB
# May 2020

library(tidyverse)
library(vegan)
library(patchwork)
source('Scripts/functions.R')

reads <- read_csv("Data/SequencedData/Plants/ProcessedData/trnL_reads.csv")
totals <- read_csv("Data/SequencedData/Plants/ProcessedData/trnL_totals.csv")
samples <- read_csv("Data/CollectionData/fecal_sample_collection.csv")

# colorblind friendly palette
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# DATA PREP #

# OTUs

# 2017, 1000, 0.01
data_2017_1000_0.01 <- filter_reads_data_trnL(samples, reads, totals, 
                                              reads_min = 1000, yr = 2017, 
                                              rel_reads_min = 0.01) %>% 
  data_prep_multivariate()
data_2017_1000_0.01[[1]] <- binarize(data_2017_1000_0.01[[1]])

dist_trnL <- metaMDS(data_2017_1000_0.01[[1]], distance = "euclidean", trymax = 250, k = 3)
(plot_2017_1000_0.01 <- NMDS_plotting_prep(data_2017_1000_0.01, dist_trnL) %>% 
    plot_NMDS_ggplot2())
(plot_2017_1000_0.01 <- plot_2017_1000_0.01 +
    ylab("Total Reads > 1000") +
    theme(axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    ggtitle("RRA > 0.01"))

# 2017, 1000, 0.05
data_2017_1000_0.05 <- filter_reads_data_trnL(samples, reads, totals, 
                                              reads_min = 1000, yr = 2017, 
                                              rel_reads_min = 0.05) %>% 
  data_prep_multivariate()
data_2017_1000_0.05[[1]] <- binarize(data_2017_1000_0.05[[1]])

dist_trnL <- metaMDS(data_2017_1000_0.05[[1]], distance = "euclidean", trymax = 250, k = 3)
(plot_2017_1000_0.05 <- NMDS_plotting_prep(data_2017_1000_0.05, dist_trnL) %>% 
    plot_NMDS_ggplot2())
(plot_2017_1000_0.05 <- plot_2017_1000_0.05 +
    theme(axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    ggtitle("RRA > 0.05"))

# 2017, 1000, 0.1
data_2017_1000_0.1 <- filter_reads_data_trnL(samples, reads, totals, 
                                             reads_min = 1000, yr = 2017, 
                                             rel_reads_min = 0.1) %>% 
  data_prep_multivariate()
data_2017_1000_0.1[[1]] <- binarize(data_2017_1000_0.1[[1]])

dist_trnL <- metaMDS(data_2017_1000_0.1[[1]], distance = "euclidean", trymax = 250, k = 3)
(plot_2017_1000_0.1 <- NMDS_plotting_prep(data_2017_1000_0.1, dist_trnL) %>% 
    plot_NMDS_ggplot2())
(plot_2017_1000_0.1 <- plot_2017_1000_0.1 +
    ggtitle("RRA > 0.1") +
    theme(axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5)))

# 2017, 2000, 0.01
data_2017_2000_0.01 <- filter_reads_data_trnL(samples, reads, totals, 
                                              reads_min = 2000, yr = 2017, 
                                              rel_reads_min = 0.01) %>% 
  data_prep_multivariate()
data_2017_2000_0.01[[1]] <- binarize(data_2017_2000_0.01[[1]])

dist_trnL <- metaMDS(data_2017_2000_0.01[[1]], distance = "euclidean", trymax = 250, k = 3)
(plot_2017_2000_0.01 <- NMDS_plotting_prep(data_2017_2000_0.01, dist_trnL) %>% 
    plot_NMDS_ggplot2())
(plot_2017_2000_0.01 <- plot_2017_2000_0.01 +
    ylab("Total Reads > 2000") +
    theme(axis.title.x = element_blank()))

# 2017, 2000, 0.05
data_2017_2000_0.05 <- filter_reads_data_trnL(samples, reads, totals, 
                                              reads_min = 2000, yr = 2017, 
                                              rel_reads_min = 0.05) %>% 
  data_prep_multivariate()
data_2017_2000_0.05[[1]] <- binarize(data_2017_2000_0.05[[1]])

dist_trnL <- metaMDS(data_2017_2000_0.05[[1]], distance = "euclidean", trymax = 250, k = 3)
(plot_2017_2000_0.05 <- NMDS_plotting_prep(data_2017_2000_0.05, dist_trnL) %>% 
    plot_NMDS_ggplot2())
(plot_2017_2000_0.05 <- plot_2017_2000_0.05 +
    theme(axis.title = element_blank()))

# 2017, 2000, 0.1
data_2017_2000_0.1 <- filter_reads_data_trnL(samples, reads, totals, 
                                              reads_min = 2000, yr = 2017, 
                                              rel_reads_min = 0.1) %>% 
  data_prep_multivariate()
data_2017_2000_0.1[[1]] <- binarize(data_2017_2000_0.1[[1]])

dist_trnL <- metaMDS(data_2017_2000_0.1[[1]], distance = "euclidean", trymax = 250, k = 3)
(plot_2017_2000_0.1 <- NMDS_plotting_prep(data_2017_2000_0.1, dist_trnL) %>% 
    plot_NMDS_ggplot2())
(plot_2017_2000_0.1 <- plot_2017_2000_0.1 +
    theme(axis.title = element_blank()))

# 2017, 5000, 0.01
data_2017_5000_0.01 <- filter_reads_data_trnL(samples, reads, totals, 
                                              reads_min = 5000, yr = 2017, 
                                              rel_reads_min = 0.01) %>% 
  data_prep_multivariate()
data_2017_5000_0.01[[1]] <- binarize(data_2017_5000_0.01[[1]])

dist_trnL <- metaMDS(data_2017_5000_0.01[[1]], distance = "euclidean", trymax = 250, k = 3)
(plot_2017_5000_0.01 <- NMDS_plotting_prep(data_2017_5000_0.01, dist_trnL) %>% 
    plot_NMDS_ggplot2())
(plot_2017_5000_0.01 <- plot_2017_5000_0.01 +
    ylab("Total Reads > 5000"))

# 2017, 5000, 0.05
data_2017_5000_0.05 <- filter_reads_data_trnL(samples, reads, totals, 
                                              reads_min = 5000, yr = 2017, 
                                              rel_reads_min = 0.05) %>% 
  data_prep_multivariate()
data_2017_5000_0.05[[1]] <- binarize(data_2017_5000_0.05[[1]])

dist_trnL <- metaMDS(data_2017_5000_0.05[[1]], distance = "euclidean", trymax = 250, k = 3)
(plot_2017_5000_0.05 <- NMDS_plotting_prep(data_2017_5000_0.05, dist_trnL) %>% 
    plot_NMDS_ggplot2())
(plot_2017_5000_0.05 <- plot_2017_5000_0.05 +
    theme(axis.title = element_blank()))

# 2017, 5000, 0.1
data_2017_5000_0.1 <- filter_reads_data_trnL(samples, reads, totals, 
                                             reads_min = 5000, yr = 2017, 
                                             rel_reads_min = 0.1) %>% 
  data_prep_multivariate()
data_2017_5000_0.1[[1]] <- binarize(data_2017_5000_0.1[[1]])

dist_trnL <- metaMDS(data_2017_5000_0.1[[1]], distance = "euclidean", trymax = 250, k = 3)
plot_2017_5000_0.1 <- NMDS_plotting_prep(data_2017_5000_0.1, dist_trnL) %>% 
    plot_NMDS_ggplot2()
(plot_2017_5000_0.1 <- plot_2017_5000_0.1 +
  theme(axis.title = element_blank()))


# put it together
patchwork <- (plot_2017_1000_0.01 + plot_2017_1000_0.05 + plot_2017_1000_0.1)/
  (plot_2017_2000_0.01 + plot_2017_2000_0.05 + plot_2017_2000_0.1)/
  (plot_2017_5000_0.01 + plot_2017_5000_0.05 + plot_2017_5000_0.1)
patchwork


# WORKING AREA ================================================================#

# remove outlier (S013067)
# data_2017_2000_0.05[[1]] <- data_2017_2000_0.05$reads_spread[-99,]
# data_2017_2000_0.05[[2]] <- data_2017_2000_0.05$sampleID[-99]
# data_2017_2000_0.05[[3]] <- data_2017_2000_0.05$groups[-99,]
