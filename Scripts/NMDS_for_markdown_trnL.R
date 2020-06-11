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


# Plot 1 #
# OTUs, 2017 ------------------------------------------------------------------#

# 2017, 1000, 0.01
dat1 <- prep_2017_allsp(samples, reads, totals, 1000, 2017, 0.01)

# 2017, 1000, 0.05
dat2 <- prep_2017_allsp(samples, reads, totals, 1000, 2017, 0.05)

# 2017, 1000, 0.005
dat3 <- prep_2017_allsp(samples, reads, totals, 1000, 2017, 0.005)

# 2017, 1000, 0.001
dat4 <- prep_2017_allsp(samples, reads, totals, 1000, 2017, 0.001)

# 2017, 2000, 0.01
dat5 <- prep_2017_allsp(samples, reads, totals, 2000, 2017, 0.01)

# 2017, 2000, 0.05
dat6 <- prep_2017_allsp(samples, reads, totals, 2000, 2017, 0.05)

# 2017, 2000, 0.005
dat7 <- prep_2017_allsp(samples, reads, totals, 2000, 2017, 0.005)

# 2017, 2000, 0.001
dat8 <- prep_2017_allsp(samples, reads, totals, 2000, 2017, 0.001)

# 2017, 5000, 0.01
dat9 <- prep_2017_allsp(samples, reads, totals, 5000, 2017, 0.01)

# 2017, 5000, 0.05
dat10 <- prep_2017_allsp(samples, reads, totals, 5000, 2017, 0.05)

# 2017, 5000, 0.005
dat11 <- prep_2017_allsp(samples, reads, totals, 5000, 2017, 0.005)

# 2017, 5000, 0.001
dat12 <- prep_2017_allsp(samples, reads, totals, 5000, 2017, 0.001)

df<- bind_rows(dat1, dat2, dat3, dat4, dat5, dat6,
               dat7, dat8, dat9, dat10, dat11, dat12)

plot1 <- ggplot(data = df[df$df == 'NMDS',], aes(x = MDS1, y = MDS2)) + 
  geom_point(aes(color = group), size = 0.5) +
  geom_path(data = df[df$df == "df_ell",], aes(x = MDS1, y = MDS2, colour = group), 
            size = 0.5) +
  facet_wrap(min_total ~ min_rel_abund, scales = "free") +
  geom_text(data = df[df$df == "NMDS.mean" & df$group == "Krat",], 
            aes(x = MDS1, y = MDS2,  
                label = .data$group[1], 
                color = .data$group[1]),
            size = 1) +
  geom_text(data = df[df$df == "NMDS.mean" & df$group == "PP_control",], 
            aes(x = MDS1, y = MDS2,  
                label = .data$group[1], 
                color = .data$group[1]),
            size = 1) +
  geom_text(data = df[df$df == "NMDS.mean" & df$group == "PP_exclosure",],
            aes(x = MDS1, y = MDS2,  
                label = .data$group[1], 
                color = .data$group[1]),
            size = 1) +
  geom_text(data = df %>% select(F.model, pval, min_total, min_rel_abund) %>% distinct(), 
            aes(x = Inf, y = Inf,
                label = paste("F.model = ", round(.data$F.model, 2),
                              "\n p = ", round(.data$pval, 4))),
            hjust = 1.1, vjust= 1.2, size = 2) +
  scale_color_manual(values = cbPalette) +
  theme_bw() +
  theme(legend.position = 'bottom',
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

#ggsave(plot1, "Plots/plot_2017_allsp_totalreads_relabund.png")

# PLOT 2 #
# OTUs, 2016 ------------------------------------------------------------------#

# 2016, 1000, 0.01
data_2016_1000_0.01 <- filter_reads_data_trnL(samples, reads, totals, 
                                              reads_min = 1000, yr = 2016, 
                                              rel_reads_min = 0.01) %>% 
  data_prep_multivariate()
data_2016_1000_0.01[[1]] <- binarize(data_2016_1000_0.01[[1]])

dist_trnL <- metaMDS(data_2016_1000_0.01[[1]], distance = "euclidean", trymax = 250, k = 3)
(plot_2016_1000_0.01 <- NMDS_plotting_prep(data_2016_1000_0.01, dist_trnL) %>% 
    plot_NMDS_ggplot2())
(plot_2016_1000_0.01 <- plot_2016_1000_0.01 +
    ylab("Total Reads > 1000") +
    theme(axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    ggtitle("RRA > 0.01"))

# 2016, 1000, 0.05
data_2016_1000_0.05 <- filter_reads_data_trnL(samples, reads, totals, 
                                              reads_min = 1000, yr = 2016, 
                                              rel_reads_min = 0.05) %>% 
  data_prep_multivariate()
data_2016_1000_0.05[[1]] <- binarize(data_2016_1000_0.05[[1]])

dist_trnL <- metaMDS(data_2016_1000_0.05[[1]], distance = "euclidean", trymax = 250, k = 3)
(plot_2016_1000_0.05 <- NMDS_plotting_prep(data_2016_1000_0.05, dist_trnL) %>% 
    plot_NMDS_ggplot2())
(plot_2016_1000_0.05 <- plot_2016_1000_0.05 +
    theme(axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    ggtitle("RRA > 0.05"))

# 2016, 1000, 0.1
data_2016_1000_0.1 <- filter_reads_data_trnL(samples, reads, totals, 
                                             reads_min = 1000, yr = 2016, 
                                             rel_reads_min = 0.1) %>% 
  data_prep_multivariate()
data_2016_1000_0.1[[1]] <- binarize(data_2016_1000_0.1[[1]])

dist_trnL <- metaMDS(data_2016_1000_0.1[[1]], distance = "euclidean", trymax = 250, k = 3)
(plot_2016_1000_0.1 <- NMDS_plotting_prep(data_2016_1000_0.1, dist_trnL) %>% 
    plot_NMDS_ggplot2())
(plot_2016_1000_0.1 <- plot_2016_1000_0.1 +
    ggtitle("RRA > 0.1") +
    theme(axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5)))

# 2016, 2000, 0.01
data_2016_2000_0.01 <- filter_reads_data_trnL(samples, reads, totals, 
                                              reads_min = 2000, yr = 2016, 
                                              rel_reads_min = 0.01) %>% 
  data_prep_multivariate()
data_2016_2000_0.01[[1]] <- binarize(data_2016_2000_0.01[[1]])

dist_trnL <- metaMDS(data_2016_2000_0.01[[1]], distance = "euclidean", trymax = 250, k = 3)
(plot_2016_2000_0.01 <- NMDS_plotting_prep(data_2016_2000_0.01, dist_trnL) %>% 
    plot_NMDS_ggplot2())
(plot_2016_2000_0.01 <- plot_2016_2000_0.01 +
    ylab("Total Reads > 2000") +
    theme(axis.title.x = element_blank()))

# 2016, 2000, 0.05
data_2016_2000_0.05 <- filter_reads_data_trnL(samples, reads, totals, 
                                              reads_min = 2000, yr = 2016, 
                                              rel_reads_min = 0.05) %>% 
  data_prep_multivariate()
data_2016_2000_0.05[[1]] <- binarize(data_2016_2000_0.05[[1]])

dist_trnL <- metaMDS(data_2016_2000_0.05[[1]], distance = "euclidean", trymax = 250, k = 3)
(plot_2016_2000_0.05 <- NMDS_plotting_prep(data_2016_2000_0.05, dist_trnL) %>% 
    plot_NMDS_ggplot2())
(plot_2016_2000_0.05 <- plot_2016_2000_0.05 +
    theme(axis.title = element_blank()))

# 2016, 2000, 0.1
data_2016_2000_0.1 <- filter_reads_data_trnL(samples, reads, totals, 
                                             reads_min = 2000, yr = 2016, 
                                             rel_reads_min = 0.1) %>% 
  data_prep_multivariate()
data_2016_2000_0.1[[1]] <- binarize(data_2016_2000_0.1[[1]])

dist_trnL <- metaMDS(data_2016_2000_0.1[[1]], distance = "euclidean", trymax = 250, k = 3)
(plot_2016_2000_0.1 <- NMDS_plotting_prep(data_2016_2000_0.1, dist_trnL) %>% 
    plot_NMDS_ggplot2())
(plot_2016_2000_0.1 <- plot_2016_2000_0.1 +
    theme(axis.title = element_blank()))

# 2016, 5000, 0.01
data_2016_5000_0.01 <- filter_reads_data_trnL(samples, reads, totals, 
                                              reads_min = 5000, yr = 2016, 
                                              rel_reads_min = 0.01) %>% 
  data_prep_multivariate()
data_2016_5000_0.01[[1]] <- binarize(data_2016_5000_0.01[[1]])

dist_trnL <- metaMDS(data_2016_5000_0.01[[1]], distance = "bray", trymax = 250, k = 3)
(plot_2016_5000_0.01 <- NMDS_plotting_prep(data_2016_5000_0.01, dist_trnL) %>% 
    plot_NMDS_ggplot2())
(plot_2016_5000_0.01 <- plot_2016_5000_0.01 +
    ylab("Total Reads > 5000"))

# 2016, 5000, 0.05
data_2016_5000_0.05 <- filter_reads_data_trnL(samples, reads, totals, 
                                              reads_min = 5000, yr = 2016, 
                                              rel_reads_min = 0.05) %>% 
  data_prep_multivariate()
data_2016_5000_0.05[[1]] <- binarize(data_2016_5000_0.05[[1]])

dist_trnL <- metaMDS(data_2016_5000_0.05[[1]], distance = "euclidean", trymax = 250, k = 3)
(plot_2016_5000_0.05 <- NMDS_plotting_prep(data_2016_5000_0.05, dist_trnL) %>% 
    plot_NMDS_ggplot2())
(plot_2016_5000_0.05 <- plot_2016_5000_0.05 +
    theme(axis.title = element_blank()))

# 2016, 5000, 0.1
data_2016_5000_0.1 <- filter_reads_data_trnL(samples, reads, totals, 
                                             reads_min = 5000, yr = 2016, 
                                             rel_reads_min = 0.1) %>% 
  data_prep_multivariate()
data_2016_5000_0.1[[1]] <- binarize(data_2016_5000_0.1[[1]])

dist_trnL <- metaMDS(data_2016_5000_0.1[[1]], distance = "euclidean", trymax = 250, k = 3)
plot_2016_5000_0.1 <- NMDS_plotting_prep(data_2016_5000_0.1, dist_trnL) %>% 
  plot_NMDS_ggplot2()
(plot_2016_5000_0.1 <- plot_2016_5000_0.1 +
    theme(axis.title = element_blank()))

# put it together
patchwork_2016 <- (plot_2016_1000_0.01 + plot_2016_1000_0.05 + plot_2016_1000_0.1)/
  (plot_2016_2000_0.01 + plot_2016_2000_0.05 + plot_2016_2000_0.1)/
  (plot_2016_5000_0.01 + plot_2016_5000_0.05 + plot_2016_5000_0.1)
patchwork_2016

# NOT READY! If using Bray-Curtis, there is an outlier! 
#   - also need to add in PP to the names
# PLOT 3 #
# OTUs, 2017, PPs only ------------------------------------------------------------------#

samples_PP <- samples %>% filter(species == 'PP')

# 2017, 1000, 0.01
data_2017_1000_0.01_PP <- filter_reads_data_trnL(samples_PP, reads, totals, 
                                              reads_min = 1000, yr = 2017, 
                                              rel_reads_min = 0.01) %>% 
  data_prep_multivariate()

dist_trnL <- metaMDS(data_2017_1000_0.01_PP[[1]], distance = "bray", trymax = 250, k = 3)
(plot_2017_1000_0.01_PP <- NMDS_plotting_prep(data_2017_1000_0.01_PP, dist_trnL) %>% 
    plot_NMDS_ggplot2())
(plot_2017_1000_0.01_PP <- plot_2017_1000_0.01_PP +
    ylab("Total Reads > 1000") +
    theme(axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    ggtitle("RRA > 0.01"))

# 2017, 1000, 0.05
data_2017_1000_0.05_PP <- filter_reads_data_trnL(samples_PP, reads, totals, 
                                              reads_min = 1000, yr = 2017, 
                                              rel_reads_min = 0.05) %>% 
  data_prep_multivariate()

dist_trnL <- metaMDS(data_2017_1000_0.05_PP[[1]], distance = "euclidean", trymax = 250, k = 3)
(plot_2017_1000_0.05_PP <- NMDS_plotting_prep(data_2017_1000_0.05_PP, dist_trnL) %>% 
    plot_NMDS_ggplot2())
(plot_2017_1000_0.05_PP <- plot_2017_1000_0.05_PP +
    theme(axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    ggtitle("RRA > 0.05"))

# 2017, 1000, 0.1
data_2017_1000_0.1_PP <- filter_reads_data_trnL(samples_PP, reads, totals, 
                                             reads_min = 1000, yr = 2017, 
                                             rel_reads_min = 0.1) %>% 
  data_prep_multivariate()

dist_trnL <- metaMDS(data_2017_1000_0.1_PP[[1]], distance = "euclidean", trymax = 250, k = 3)
(plot_2017_1000_0.1_PP <- NMDS_plotting_prep(data_2017_1000_0.1_PP, dist_trnL) %>% 
    plot_NMDS_ggplot2())
(plot_2017_1000_0.1_PP <- plot_2017_1000_0.1_PP +
    ggtitle("RRA > 0.1") +
    theme(axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5)))

# 2017, 2000, 0.01
data_2017_2000_0.01_PP <- filter_reads_data_trnL(samples_PP, reads, totals, 
                                              reads_min = 2000, yr = 2017, 
                                              rel_reads_min = 0.01) %>% 
  data_prep_multivariate()

dist_trnL <- metaMDS(data_2017_2000_0.01_PP[[1]], distance = "euclidean", trymax = 250, k = 3)
(plot_2017_2000_0.01_PP <- NMDS_plotting_prep(data_2017_2000_0.01_PP, dist_trnL) %>% 
    plot_NMDS_ggplot2())
(plot_2017_2000_0.01_PP <- plot_2017_2000_0.01_PP +
    ylab("Total Reads > 2000") +
    theme(axis.title.x = element_blank()))

# 2017, 2000, 0.05
data_2017_2000_0.05_PP <- filter_reads_data_trnL(samples_PP, reads, totals, 
                                              reads_min = 2000, yr = 2017, 
                                              rel_reads_min = 0.05) %>% 
  data_prep_multivariate()

dist_trnL <- metaMDS(data_2017_2000_0.05_PP[[1]], distance = "euclidean", trymax = 250, k = 3)
(plot_2017_2000_0.05_PP <- NMDS_plotting_prep(data_2017_2000_0.05_PP, dist_trnL) %>% 
    plot_NMDS_ggplot2())
(plot_2017_2000_0.05_PP <- plot_2017_2000_0.05_PP +
    theme(axis.title = element_blank()))

# 2017, 2000, 0.1
data_2017_2000_0.1_PP <- filter_reads_data_trnL(samples_PP, reads, totals, 
                                             reads_min = 2000, yr = 2017, 
                                             rel_reads_min = 0.1) %>% 
  data_prep_multivariate()

dist_trnL <- metaMDS(data_2017_2000_0.1_PP[[1]], distance = "euclidean", trymax = 250, k = 3)
(plot_2017_2000_0.1_PP <- NMDS_plotting_prep(data_2017_2000_0.1_PP, dist_trnL) %>% 
    plot_NMDS_ggplot2())
(plot_2017_2000_0.1_PP <- plot_2017_2000_0.1_PP +
    theme(axis.title = element_blank()))

# 2017, 5000, 0.01
data_2017_5000_0.01_PP <- filter_reads_data_trnL(samples_PP, reads, totals, 
                                              reads_min = 5000, yr = 2017, 
                                              rel_reads_min = 0.01) %>% 
  data_prep_multivariate()

dist_trnL <- metaMDS(data_2017_5000_0.01_PP[[1]], distance = "euclidean", trymax = 250, k = 3)
(plot_2017_5000_0.01_PP <- NMDS_plotting_prep(data_2017_5000_0.01_PP, dist_trnL) %>% 
    plot_NMDS_ggplot2())
(plot_2017_5000_0.01_PP <- plot_2017_5000_0.01_PP +
    ylab("Total Reads > 5000"))

# 2017, 5000, 0.05
data_2017_5000_0.05_PP <- filter_reads_data_trnL(samples_PP, reads, totals, 
                                              reads_min = 5000, yr = 2017, 
                                              rel_reads_min = 0.05) %>% 
  data_prep_multivariate()

dist_trnL <- metaMDS(data_2017_5000_0.05_PP[[1]], distance = "euclidean", trymax = 250, k = 3)
(plot_2017_5000_0.05_PP <- NMDS_plotting_prep(data_2017_5000_0.05_PP, dist_trnL) %>% 
    plot_NMDS_ggplot2())
(plot_2017_5000_0.05_PP <- plot_2017_5000_0.05_PP +
    theme(axis.title = element_blank()))

# 2017, 5000, 0.1
data_2017_5000_0.1_PP <- filter_reads_data_trnL(samples_PP, reads, totals, 
                                             reads_min = 5000, yr = 2017, 
                                             rel_reads_min = 0.1) %>% 
  data_prep_multivariate()

dist_trnL <- metaMDS(data_2017_5000_0.1_PP[[1]], distance = "euclidean", trymax = 250, k = 3)
plot_2017_5000_0.1_PP <- NMDS_plotting_prep(data_2017_5000_0.1_PP, dist_trnL) %>% 
  plot_NMDS_ggplot2()
(plot_2017_5000_0.1_PP <- plot_2017_5000_0.1_PP +
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

data_2017_1000_0.01[[1]] <- data_2017_1000_0.01[[1]][!(row.names(data_2017_1000_0.01[[1]]) %in% c("S010049")),]

# scree plot
goeveg::dimcheckMDS(data_2017_2000_0.05[[1]], distance = "bray", k = 6, trymax = 50)


prep1 <- NMDS_plotting_prep(data_2017_1000_0.01, dist_trnL)
prep2 <- NMDS_plotting_prep(data_2017_1000_0.005, dist_trnL)


prep1[[1]]$df <- "NMDS"
prep1[[2]]$df <- "NMDS.mean"
prep1[[3]]$df <- "df_ell"
prep1[[3]] <- prep1[[3]] %>% 
  rename("MDS1" = NMDS1, "MDS2" = NMDS2)
df1 <- bind_rows(prep1[[1]], prep1[[2]], prep1[[3]])
df1$F.model <- prep1[[4]]$aov.tab$F.Model[1]
df1$pval <- prep1[[4]]$aov.tab$`Pr(>F)`[1]
df1$min_total <- 1000
df1$min_rel_abund <- 0.01

prep2[[1]]$df <- "NMDS"
prep2[[2]]$df <- "NMDS.mean"
prep2[[3]]$df <- "df_ell"
prep2[[3]] <- prep2[[3]] %>% 
  rename("MDS1" = NMDS1, "MDS2" = NMDS2)
df2 <- bind_rows(prep2[[1]], prep2[[2]], prep2[[3]])
df2$F.model <- prep2[[4]]$aov.tab$F.Model[1]
df2$pval <- prep2[[4]]$aov.tab$`Pr(>F)`[1]
df2$min_total <- 1000
df2$min_rel_abund <- 0.005

df <- bind_rows(df1, df2)

ggplot(data = df[df$df == 'NMDS',], aes(x = MDS1, y = MDS2)) + 
  geom_point(aes(color = group), size = 0.5) +
  geom_path(data = df[df$df == "df_ell",], aes(x = MDS1, y = MDS2, colour = group), 
            size = 0.5) +
  facet_grid(~ min_rel_abund) +
  geom_text(data = df[df$df == "NMDS.mean" & df$group == "Krat",], 
            aes(x = MDS1, y = MDS2,  
                label = .data$group[1], 
                color = .data$group[1]),
            size = 1) +
  geom_text(data = df[df$df == "NMDS.mean" & df$group == "PP_control",], 
            aes(x = MDS1, y = MDS2,  
                label = .data$group[1], 
                color = .data$group[1]),
            size = 1) +
  geom_text(data = df[df$df == "NMDS.mean" & df$group == "PP_exclosure",],
            aes(x = MDS1, y = MDS2,  
                label = .data$group[1], 
                color = .data$group[1]),
            size = 1) +
  geom_text(data = df %>% select(F.model, pval, min_rel_abund) %>% distinct(), 
            aes(x = Inf, y = Inf,
                label = paste("F.model = ", round(.data$F.model, 2),
                "\n p = ", round(.data$pval, 4))),
            hjust = 1.1, vjust= 1.2, size = 2) +
  scale_color_manual(values = cbPalette) +
  theme_bw() +
  theme(legend.position = 'none',
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 


df %>% 
  select(F.model, pval, min_total, min_rel_abund) %>% 
  distinct()

  # annotate(geom = "text", x = Inf, y = Inf, hjust = 1.1, vjust= 1.2, 
  #          size = 2, group = .data$min_rel_abund,
  #          label = paste("atop(' F.model = '*", round(df$F.model[1], 2),"
  #                        ,' p = '*", round(df$pval[1], 4),")"), parse=T)
