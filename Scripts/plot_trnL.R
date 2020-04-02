# First Pass at Plotting 
# Fecal samples, trnL
# March 2020 EKB

# LIBRARIES, SOURCE CODE & DATA #

library(tidyverse)
library(vegan)
source("Scripts/functions.R")

reads <- read_csv("Data/SequencedData/Plants/ProcessedData/trnL_reads.csv")
totals <- read_csv("Data/SequencedData/Plants/ProcessedData/trnL_totals.csv")
samples <- read_csv("Data/CollectionData/fecal_sample_collection.csv")

# colorblind friendly palette
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# DATA PREP #

test_list_noyear <- filter_reads_data(samples, reads, totals)
test_list_2017 <- filter_reads_data(samples, reads, totals, yr = 2017, rel_reads_min = 0.005)
test_list_2016 <- filter_reads_data(samples, reads, totals, yr = 2016)
test_list_noyear_0.005 <- filter_reads_data(samples, reads, totals, rel_reads_min = 0.005)

test_list2 <- data_prep_for_NMDS(test_list_2017)

# NMDS ANALYSIS #

# check coefficient of variation to see if standardization is needed
# cSums <- colSums(reads_spread)
# Sdev <- sd(cSums)
# M <- mean(cSums)
# CV <- (Sdev/M)*100
# 
# # scale the data (Z-standardization)
# reads_scaled <- scale(reads_spread)
# 
# # detect outliers
# out<-function(x){
#   lier<-x[abs(x)>3]
#   return(lier)
# }
# 
# apply(reads_scaled, 2, out)


# PLOTTING #

#### need to pull out dataframes from test_list2 (or equivalent) 
#    and set "reads_spread" and "groups"



# scree plot
goeveg::dimcheckMDS(test_list2[[1]], distance = "euclidean", k = 6, trymax = 50)

# use vegan package to run NMDS
dist_trnL <- metaMDS(test_list2[[1]], distance = "euclidean", trymax = 250, k = 3, 
                     noshare = 0.2)
dist_trnL <- metaMDS(test_list2[[1]], distance = "euclidean", trymax = 50, k = 3, 
                     noshare = 0.2, previous.best = dist_trnL)
stressplot(dist_trnL)

# make dataframe with with MDS values and grouping variable
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

# run PERMANOVA
group = as.matrix(groups$group)
permtrnL <- adonis(reads_spread ~ group, permutations = 10000)
permtrnL
hist(permtrnL$f.perms)
points(permtrnL$aov.tab$F.Model[1], 0 , pch = 19, col = "red", bg = "red", cex = 2)
