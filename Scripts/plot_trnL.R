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

# add plot type to fecal collection data
# add group for plotting
# and remove samples that were part of the trap/bait test
samples <- add_plot_type(samples) %>% 
  add_plotting_group() %>% 
  filter(is.na(notes), year == 2017)

# select only fecal samples
fecal_id <- samples$vial_barcode

# add totals to reads df, filter out fecal samples and small totals
reads <- full_join(reads, totals)
reads <- reads %>% 
  filter(SampleID %in% fecal_id) %>% 
  filter(Total_Reads > 2000) %>% 
  mutate(Rel_Reads = Reads/Total_Reads)

reads <- select(reads, SampleID, OTU, Rel_Reads)
reads_spread <- pivot_wider(reads, names_from = OTU, values_from = Rel_Reads)
reads_spread[is.na(reads_spread)] = 0
reads_spread <- reads_spread %>% tibble::column_to_rownames("SampleID")

sampleID <- intersect(reads$SampleID, samples$vial_barcode)
groups <- samples %>% 
  filter(vial_barcode %in% sampleID)

# NMDS ANALYSIS #

# check coefficient of variation to see if standardization is needed
cSums <- colSums(reads_spread)
Sdev <- sd(cSums)
M <- mean(cSums)
CV <- (Sdev/M)*100

# scale the data (Z-standardization)
reads_scaled <- scale(reads_spread)

# detect outliers
out<-function(x){
  lier<-x[abs(x)>3]
  return(lier)
}

apply(reads_scaled, 2, out)


# PLOTTING #

# scree plot
goeveg::dimcheckMDS(reads_scaled, distance = "euclidean", k = 6, trymax = 50)

# use vegan package to run NMDS
dist_trnL <- metaMDS(reads_spread, distance = "euclidean", trymax = 50, k = 3, 
                     noshare = TRUE)
dist_trnL <- metaMDS(reads_spread, distance = "euclidean", trymax = 50, k = 3, 
                     noshare = TRUE, previous.best = dist_trnL)
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
