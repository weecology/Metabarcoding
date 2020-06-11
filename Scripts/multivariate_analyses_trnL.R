# Multivariate Attempts
# EKB
# May 2020

# LIBRARIES AND DATA #

library(tidyverse)
library(vegan)
source("Scripts/functions.R")

reads <- read_csv("Data/SequencedData/Plants/ProcessedData/trnL_reads.csv")
totals <- read_csv("Data/SequencedData/Plants/ProcessedData/trnL_totals.csv")
samples <- read_csv("Data/CollectionData/fecal_sample_collection.csv")

# colorblind friendly palette
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# DATA PREP #

test_list_2017 <- filter_reads_data_trnL(samples, reads, totals, yr = 2017, rel_reads_min = 0.01)
test_list_2016 <- filter_reads_data_trnL(samples, reads, totals, yr = 2016)
test_list_noyear_0.005 <- filter_reads_data_trnL(samples, reads, totals, rel_reads_min = 0.005)

# check this to see if it works
test_list_2017[[1]] <- test_list2[[1]] %>% select(-OTU150)

test_list2 <- data_prep_multivariate(test_list_2017)

# remove outlier (S013067)
test_list2[[1]] <- test_list2$reads_spread[-105,]
test_list2[[2]] <- test_list2$sampleID[-105]
test_list2[[3]] <- test_list2$groups[-105,]

test_list2_binary <- binarize(test_list2[[1]])

#===============================================================================
# Non-metric Dimensional Scaling (NMDS) #
#===============================================================================

# scree plot
goeveg::dimcheckMDS(test_list2[[1]], distance = "bray", k = 6, trymax = 50)

# use vegan package to run NMDS
dist_trnL <- metaMDS(test_list2[[1]], distance = "bray", trymax = 250, k = 3)
dist_trnL_converge <- run_metaMDS_til_converge(test_list2[[1]], dist_trnL, "euclidean", 3)

stressplot(dist_trnL)

# make dataframe with with MDS values and grouping variable
groups <- test_list2[[3]]

NMDS <- data.frame(MDS1 = dist_trnL$points[,1], 
                   MDS2 = dist_trnL$points[,2], 
                   group = groups$group)

# get mean point for each group
NMDS.mean <- aggregate(NMDS[,1:2], list(group = groups$group), mean)

# save results of ordiellipse() as an object
plot(dist_trnL$points)
ord <- ordiellipse(dist_trnL, 
                   groups$group, 
                   display = "sites", 
                   kind = "se", 
                   conf = 0.95, 
                   label = T)

# plot using ggplot 
df_ell <- data.frame()
for(g in levels(NMDS$group)) {
  df_ell <-
    rbind(df_ell, cbind(as.data.frame(with(
      NMDS[NMDS$group == g,],
      vegan:::veganCovEllipse(ord[[g]]$cov, ord[[g]]$center, ord[[g]]$scale)
    )),
    group = g))
}

ggplot(data = NMDS, aes(x = MDS1, y = MDS2)) + 
  geom_point(aes(color = group)) +
  geom_path(data = df_ell, aes(x = NMDS1, y = NMDS2, colour = group), size = 1) +
  geom_text(aes(x = NMDS.mean$MDS1[1], y = NMDS.mean$MDS2[1], label = NMDS.mean$group[1], color = NMDS.mean$group[1])) +
  geom_text(aes(x = NMDS.mean$MDS1[2], y = NMDS.mean$MDS2[2], label = NMDS.mean$group[2], color = NMDS.mean$group[2])) +
  geom_text(aes(x = NMDS.mean$MDS1[3], y = NMDS.mean$MDS2[3], label = NMDS.mean$group[3], color = NMDS.mean$group[3])) +
  scale_color_manual(values = cbPalette) +
  theme_bw() +
  theme(legend.position = 'none')

# NMDS NOTES

# 2017 is the one providing issues
# euclidean (maybe bray-curtis) for rel abundance; jaccard or bray for binary
# at rel_read_abund = 0.01, k = 3 on binary 2017 data, we get convergence
#   - with outlier (S013067) removed
#   - it shows basically the exact same thing as rel abundance but with convergence!







#===============================================================================
# Principle Components Analysis (PCA)
#===============================================================================




#===============================================================================
# PCoA
#===============================================================================

#===============================================================================
# Detrended Correspondance Analysis (DCA)
#===============================================================================

