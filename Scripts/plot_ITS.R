# First Pass at Plotting 
# Fecal samples, ITS
# April 2020 EKB

# LIBRARIES, SOURCE CODE & DATA #

library(tidyverse)
library(vegan)
source("Scripts/functions.R")

reads <- read_csv("Data/SequencedData/Plants/ProcessedData/ITS2_reads.csv")
totals <- read_csv("Data/SequencedData/Plants/ProcessedData/ITS2_totals.csv")
samples <- read_csv("Data/CollectionData/fecal_sample_collection.csv")

# colorblind friendly palette
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# NOTE!
#   - these data come from multiple datasets, so the OTUs don't fully match up
#   - that has NOT yet been dealt with, so you need to figure that out

# DATA PREP #

test_list2 <- filter_reads_data(samples, reads, totals, 
                                yr = 2017, 
                                reads_min = 2000,
                                rel_reads_min = 0.005) %>% 
  data_prep_for_NMDS()

# NMDS ANALYSIS #

# scree plot
goeveg::dimcheckMDS(test_list2[[1]], distance = "bray", k = 6, trymax = 50)

# use vegan package to run NMDS
dist_ITS <- metaMDS(test_list2[[1]], distance = "bray", trymax = 250, k = 3, 
                     noshare = 0.2)
dist_ITS <- metaMDS(test_list2[[1]], distance = "bray", trymax = 250, k = 3, 
                     noshare = 0.2, previous.best = dist_ITS)
stressplot(dist_ITS)

# make dataframe with with MDS values and grouping variable
NMDS <- data.frame(MDS1 = dist_ITS$points[,1], 
                   MDS2 = dist_ITS$points[,2], 
                   group = test_list2[[3]]$group)

# get mean point for each group
NMDS.mean <- aggregate(NMDS[,1:2], list(group = test_list2[[3]]$group), mean)

# save results of ordiellipse() as an object
plot(dist_ITS$points)
ord <- ordiellipse(dist_ITS, 
                   test_list2[[3]]$group, 
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
