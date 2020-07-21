# NMDS trnL Markdown Plotting: WeeTU
# EKB
# June 2020

# LIBRARIES and DATA #

library(tidyverse)
library(vegan)
source('Scripts/functions.R')

reads <- read_csv("Data/SequencedData/Plants/ProcessedData/trnL_reads_WeeTU.csv")
totals <- read_csv("Data/SequencedData/Plants/ProcessedData/trnL_totals.csv")
samples <- read_csv("Data/CollectionData/fecal_sample_collection.csv")
OTU_WTU_key <- read_csv("Data/SequencedData/Plants/ProcessedData/trnL_OTU_WTU_key.csv")

# colorblind friendly palette
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# DATA PREP #

# only PPs
samples_PP <- samples %>% filter(species == 'PP')


# WeeTUS #======================================================================

# trnL, 2017, all species 

reads_genus <- summarize_trnL_by_WeeTU(reads, "WTU.genus", WTU.genus)
dat1 <- prep_2017_allsp_relabund_WTU_genus(samples, reads_genus, totals, OTU_WTU_key, "genus", 2000, 2017, 0.01)
dat2 <- prep_2017_allsp_relabund_WTU_genus(samples, reads_genus, totals, OTU_WTU_key, "genus", 2000, 2017, 0.05)
dat3 <- prep_2017_allsp_relabund_WTU_genus(samples, reads_genus, totals, OTU_WTU_key, "genus", 2000, 2017, 0.005)
dat4 <- prep_2017_allsp_relabund_WTU_genus(samples, reads_genus, totals, OTU_WTU_key, "genus", 2000, 2017, 0.001)

reads_family <- summarize_trnL_by_WeeTU(reads, "WTU.family", WTU.family)
dat5 <- prep_2017_allsp_relabund_WTU_family(samples, reads_family, totals, OTU_WTU_key, "family", 2000, 2017, 0.01)
dat6 <- prep_2017_allsp_relabund_WTU_family(samples, reads_family, totals, OTU_WTU_key, "family", 2000, 2017, 0.05)
dat7 <- prep_2017_allsp_relabund_WTU_family(samples, reads_family, totals, OTU_WTU_key, "family", 2000, 2017, 0.005)
dat8 <- prep_2017_allsp_relabund_WTU_family(samples, reads_family, totals, OTU_WTU_key, "family", 2000, 2017, 0.001)

df <- bind_rows(dat1, dat2, dat3, dat4, dat5, dat6,
                dat7, dat8)

(plot1 <- ggplot(data = df[df$df == 'NMDS',], aes(x = MDS1, y = MDS2)) + 
    geom_point(aes(color = group), size = 0.5) +
    geom_path(data = df[df$df == "df_ell",], aes(x = MDS1, y = MDS2, colour = group), 
              size = 0.5) +
    facet_wrap(sum_taxa ~ min_rel_abund, scales = "free", nrow = 2, ncol = 4) +
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
    geom_text(data = df %>% select(F.model, pval, sum_taxa, min_rel_abund) %>% distinct(), 
              aes(x = Inf, y = Inf,
                  label = paste("F.model = ", round(.data$F.model, 2),
                                "\n p = ", round(.data$pval, 4))),
              hjust = 1.1, vjust= 1.2, size = 2) +
    scale_color_manual(values = cbPalette) +
    ggtitle("trnL: Spring 2017, WTU") +
    theme_bw() +
    theme(legend.position = 'bottom',
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())) 

#ggsave("Plots/trnL_2017_allsp_relabund_WTU.png", plot1)

# trnL, 2017, PP only

reads_genus <- summarize_trnL_by_WeeTU(reads, "WTU.genus", WTU.genus)
dat1 <- prep_2017_PPonly_relabund_WTU_genus(samples_PP, reads_genus, totals, OTU_WTU_key, "genus", 2000, 2017, 0.01)
dat2 <- prep_2017_PPonly_relabund_WTU_genus(samples_PP, reads_genus, totals, OTU_WTU_key, "genus", 2000, 2017, 0.05)
dat3 <- prep_2017_PPonly_relabund_WTU_genus(samples_PP, reads_genus, totals, OTU_WTU_key, "genus", 2000, 2017, 0.005)
dat4 <- prep_2017_PPonly_relabund_WTU_genus(samples_PP, reads_genus, totals, OTU_WTU_key, "genus", 2000, 2017, 0.001)

reads_family <- summarize_trnL_by_WeeTU(reads, "WTU.family", WTU.family)
  # I think the summary is just too coarse here
dat5 <- prep_2017_PPonly_relabund_WTU_family(samples_PP, reads_family, totals, OTU_WTU_key, "family", 2000, 2017, 0.01)
dat6 <- prep_2017_PPonly_relabund_WTU_family(samples_PP, reads_family, totals, OTU_WTU_key, "family", 2000, 2017, 0.05)
dat7 <- prep_2017_PPonly_relabund_WTU_family(samples_PP, reads_family, totals, OTU_WTU_key, "family", 2000, 2017, 0.005)
dat8 <- prep_2017_PPonly_relabund_WTU_family(samples_PP, reads_family, totals, OTU_WTU_key, "family", 2000, 2017, 0.001)

df <- bind_rows(dat1, dat2, dat3, dat4, dat5, dat6,
                dat7, dat8)

(plot2 <- ggplot(data = df[df$df == 'NMDS',], aes(x = MDS1, y = MDS2)) + 
    geom_point(aes(color = group), size = 0.5) +
    geom_path(data = df[df$df == "df_ell",], aes(x = MDS1, y = MDS2, colour = group), 
              size = 0.5) +
    facet_wrap(sum_taxa ~ min_rel_abund, scales = "free", nrow = 2, ncol = 4) +
    # geom_text(data = df[df$df == "NMDS.mean" & df$group == "Krat",], 
    #           aes(x = MDS1, y = MDS2,  
    #               label = .data$group[1], 
    #               color = .data$group[1]),
    #           size = 1) +
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
    geom_text(data = df %>% select(F.model, pval, sum_taxa, min_rel_abund) %>% distinct(), 
              aes(x = Inf, y = Inf,
                  label = paste("F.model = ", round(.data$F.model, 2),
                                "\n p = ", round(.data$pval, 4))),
              hjust = 1.1, vjust= 1.2, size = 2) +
    scale_color_manual(values = cbPalette) +
    ggtitle("trnL: Spring 2017, WTU") +
    theme_bw() +
    theme(legend.position = 'bottom',
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())) 

#ggsave("Plots/trnL_2017_PPonly_relabund_WTU.png", plot2)

# trnL, 2016, all species 

reads_genus <- summarize_trnL_by_WeeTU(reads, "WTU.genus", WTU.genus)
dat1 <- prep_2016_allsp_relabund_WTU(samples, reads_genus, totals, OTU_WTU_key, "genus", 2000, 2016, 0.01)
dat2 <- prep_2016_allsp_relabund_WTU(samples, reads_genus, totals, OTU_WTU_key, "genus", 2000, 2016, 0.05)
dat3 <- prep_2016_allsp_relabund_WTU(samples, reads_genus, totals, OTU_WTU_key, "genus", 2000, 2016, 0.005)
dat4 <- prep_2016_allsp_relabund_WTU(samples, reads_genus, totals, OTU_WTU_key, "genus", 2000, 2016, 0.001)

reads_family <- summarize_trnL_by_WeeTU(reads, "WTU.family", WTU.family)
dat5 <- prep_2016_allsp_relabund_WTU(samples, reads_family, totals, OTU_WTU_key, "family", 2000, 2016, 0.01)
dat6 <- prep_2016_allsp_relabund_WTU(samples, reads_family, totals, OTU_WTU_key, "family", 2000, 2016, 0.05)
dat7 <- prep_2016_allsp_relabund_WTU(samples, reads_family, totals, OTU_WTU_key, "family", 2000, 2016, 0.005)
dat8 <- prep_2016_allsp_relabund_WTU(samples, reads_family, totals, OTU_WTU_key, "family", 2000, 2016, 0.001)

df <- bind_rows(dat1, dat2, dat3, dat4, dat5, dat6,
                dat7, dat8)

(plot3 <- ggplot(data = df[df$df == 'NMDS',], aes(x = MDS1, y = MDS2)) + 
    geom_point(aes(color = group), size = 0.5) +
    geom_path(data = df[df$df == "df_ell",], aes(x = MDS1, y = MDS2, colour = group), 
              size = 0.5) +
    facet_wrap(sum_taxa ~ min_rel_abund, scales = "free", nrow = 2, ncol = 4) +
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
    geom_text(data = df %>% select(F.model, pval, sum_taxa, min_rel_abund) %>% distinct(), 
              aes(x = Inf, y = Inf,
                  label = paste("F.model = ", round(.data$F.model, 2),
                                "\n p = ", round(.data$pval, 4))),
              hjust = 1.1, vjust= 1.2, size = 2) +
    scale_color_manual(values = cbPalette) +
    ggtitle("trnL: Spring 2017, WTU") +
    theme_bw() +
    theme(legend.position = 'bottom',
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())) 

#ggsave("Plots/trnL_2016_allsp_relabund_WTU.png", plot3)

# trnL, 2016, PP only

reads_genus <- summarize_trnL_by_WeeTU(reads, "WTU.genus", WTU.genus)
dat1 <- prep_2016_PPonly_relabund_WTU(samples_PP, reads_genus, totals, OTU_WTU_key, "genus", 2000, 2016, 0.01)
dat2 <- prep_2016_PPonly_relabund_WTU(samples_PP, reads_genus, totals, OTU_WTU_key, "genus", 2000, 2016, 0.05)
dat3 <- prep_2016_PPonly_relabund_WTU(samples_PP, reads_genus, totals, OTU_WTU_key, "genus", 2000, 2016, 0.005)
dat4 <- prep_2016_PPonly_relabund_WTU(samples_PP, reads_genus, totals, OTU_WTU_key, "genus", 2000, 2016, 0.001)

reads_family <- summarize_trnL_by_WeeTU(reads, "WTU.family", WTU.family)
# I think the summary is just too coarse here
dat5 <- prep_2016_PPonly_relabund_WTU(samples_PP, reads_family, totals, OTU_WTU_key, "family", 2000, 2016, 0.01)
dat6 <- prep_2016_PPonly_relabund_WTU(samples_PP, reads_family, totals, OTU_WTU_key, "family", 2000, 2016, 0.05)
dat7 <- prep_2016_PPonly_relabund_WTU(samples_PP, reads_family, totals, OTU_WTU_key, "family", 2000, 2016, 0.005)
dat8 <- prep_2016_PPonly_relabund_WTU(samples_PP, reads_family, totals, OTU_WTU_key, "family", 2000, 2016, 0.001)

df <- bind_rows(dat1, dat2, dat3, dat4, dat5, dat6,
                dat7, dat8)

(plot4 <- ggplot(data = df[df$df == 'NMDS',], aes(x = MDS1, y = MDS2)) + 
    geom_point(aes(color = group), size = 0.5) +
    geom_path(data = df[df$df == "df_ell",], aes(x = MDS1, y = MDS2, colour = group), 
              size = 0.5) +
    facet_wrap(sum_taxa ~ min_rel_abund, scales = "free", nrow = 2, ncol = 4) +
    # geom_text(data = df[df$df == "NMDS.mean" & df$group == "Krat",], 
    #           aes(x = MDS1, y = MDS2,  
    #               label = .data$group[1], 
    #               color = .data$group[1]),
    #           size = 1) +
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
    geom_text(data = df %>% select(F.model, pval, sum_taxa, min_rel_abund) %>% distinct(), 
              aes(x = Inf, y = Inf,
                  label = paste("F.model = ", round(.data$F.model, 2),
                                "\n p = ", round(.data$pval, 4))),
              hjust = 1.1, vjust= 1.2, size = 2) +
    scale_color_manual(values = cbPalette) +
    ggtitle("trnL: Spring 2017, WTU") +
    theme_bw() +
    theme(legend.position = 'bottom',
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())) 

ggsave("Plots/trnL_2016_PPonly_relabund_WTU.png", plot4)

# WORKING AREA ================================================================#

data <- filter_reads_data_WeeTU_trnL(samples_PP, 
                                     reads_genus, 
                                     totals,
                                     OTU_WTU_key,
                                     yr = 2016, 
                                     reads_min = 2000, 
                                     rel_reads_min = 0.01) %>% 
  data_prep_multivariate_WTU()
#data[[1]] <- binarize(data[[1]])

# remove outliers
# group 1: "S008810", "S010014", "S013043"
# group 2: group 1 + "S010063", "S010044", "S010012"
# data[[1]] <-
#   data[[1]][!(row.names(data[[1]]) %in% c("S008804", "S013025", "S013017",
#                                           "S010029", "S010068")),]
# data[[2]] <-
#   data[[2]][!data[[2]] %in% c("S008804", "S013025", "S013017",
#                               "S010029", "S010068")]
# data[[3]] <-
#   data[[3]][!(data[[3]]$vial_barcode) %in% c("S008804", "S013025", "S013017",
#                                              "S010029", "S010068"),]

dist_trnL <- metaMDS(data[[1]], distance = "bray", trymax = 250, k = 3)
dist_trnL$points

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

# scree plot
goeveg::dimcheckMDS(data[[1]], distance = "bray", k = 6, trymax = 50)