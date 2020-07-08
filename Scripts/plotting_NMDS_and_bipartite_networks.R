# Figure Plotting Examples
#     * just some code that will be helpful when plotting for pubs
#     * such as pairwise perMANOVAs and bipartite network graphs
# EKB 
# July 2020

# LIBRARIES and DATA #

library(tidyverse)
library(vegan)
library(EcolUtils) # pairwise perMANOVA function
source('Scripts/functions.R')

reads <- read_csv("Data/SequencedData/Plants/ProcessedData/ITS2_reads_WeeTU.csv")
totals <- read_csv("Data/SequencedData/Plants/ProcessedData/ITS2_totals.csv")
samples <- read_csv("Data/CollectionData/fecal_sample_collection.csv")

# colorblind friendly palette
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# NMDS PLOTTING #---------------------------------------------------------------

data <- filter_reads_data_ITS2(samples,
                               reads,
                               totals,
                               reads_min = 2000,
                               yr = 2016,
                               rel_reads_min = 0.01) %>%
  data_prep_multivariate()
data[[1]] <- binarize(data[[1]])

# remove outliers
data[[1]] <- 
  data[[1]][!(row.names(data[[1]]) %in% c("S008824")),]
data[[2]] <- 
  data[[2]][!data[[2]] %in% c("S008824")]
data[[3]] <- 
  data[[3]][!(data[[3]]$vial_barcode) %in% c("S008824"),]

# get distance points and matrix
dist_trnL <- metaMDS(data[[1]], distance = "bray", trymax = 250, k = 3)
dist_trnL$points

dist_matrix <- metaMDSredist(dist_trnL)

# make dataframe for plotting
groups <- data[[3]]
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

# run perMANOVA
group = as.matrix(groups$group)
perMANOVA_output <- adonis(data[[1]] ~ group, permutations = 10000)
pairwise_perMANOVA <- adonis.pair(dist.mat = dist_matrix, 
                                  Factor = as.factor(groups$group),
                                  nper = 10000)

# BIPARTITE PLOTTING #----------------------------------------------------------

library(bipartite)

data2 <- filter_reads_data_ITS2(samples,
                               reads,
                               totals,
                               reads_min = 2000,
                               yr = 2016,
                               rel_reads_min = 0.01) 

# remove outliers
data2[[1]] <- 
  data2[[1]][!(data2[[1]]$vial_barcode %in% c("S008824")),]
data2[[2]] <- 
  data2[[2]][!data2[[2]] %in% c("S008824")]
data2[[3]] <- 
  data2[[3]][!(data2[[3]]$SampleID) %in% c("S008824"),]

samples2 <- data2[[1]]
reads2 <- data2[[3]]

# convert to appropriate format for NMDS
reads2 <- select(reads2, SampleID, OTU, Rel_Reads)
reads_spread <- pivot_wider(reads2, names_from = OTU, values_from = Rel_Reads)
reads_spread[is.na(reads_spread)] = 0
reads_spread <- left_join(reads_spread, 
                          select(samples2, vial_barcode, group), 
                          by = c("SampleID" = "vial_barcode"))
reads_spread <- reads_spread %>% 
  select(group, starts_with("OTU"), -SampleID)

reads_spread[,2:48] <- binarize(reads_spread[,2:48])

reads_spread2 <- reads_spread %>% 
  group_by(group) %>% 
  summarise_each(funs(sum)) %>% 
  column_to_rownames("group")

grouplevel(reads_spread2, level = "lower") 
  # lower is actually rodents for plotting reasons

visweb(reads_spread2)
plotweb(reads_spread2,
        method = "cca", labsize = 1, 
        ybig = 1,  y.width.low = 0.1, y.width.high = 0.1, 
        col.interaction="grey80", 
        #col.high = c("mediumblue","lightblue4","salmon1", "salmon2", "salmon", "salmon3", "salmon4") ,
        col.low=c("#E69F00", "#56B4E9", "#009E73"),  
        text.rot=90, 
        #text.high.col=c("mediumblue","lightblue4","salmon1", "salmon2", "salmon", "salmon3", "salmon4"), 
        text.low.col="grey45", 
        adj.high=NULL, adj.low=NULL, plot.axes = TRUE,
        low.y=0.5, high.y=1.5, add=FALSE, y.lim=NULL, x.lim=NULL, 
        low.plot=TRUE, high.plot=TRUE, high.xoff = 0, low.xoff = 0, 
        high.lab.dis = NULL, low.lab.dis = NULL, abuns.type="additional")

# alternative package
library(bipartiteD3)

reads3 <- left_join(reads2, 
                    select(samples2, vial_barcode, group), 
                    by = c("SampleID" = "vial_barcode"))
reads3 <- reads3 %>% 
  select(-SampleID) %>% 
  group_by(group, OTU) %>% 
  summarise(freq = n()) %>% 
  ungroup()
reads3$webID <- "web1"
reads3 <- rename(reads3, 
                 higher = group, lower = OTU) %>% 
  select(higher, lower, webID, freq)

# has to be a dataframe instead of a tibble
test_web <- bipartite::frame2webs(as.data.frame(reads3))
bipartite_D3(test_web)
