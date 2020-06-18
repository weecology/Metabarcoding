# NMDS ITS2 Markdown Plotting
# EKB
# June 2020

# LIBRARIES and DATA #

library(tidyverse)
library(vegan)
source('Scripts/functions.R')

reads <- read_csv("Data/SequencedData/Plants/ProcessedData/ITS2_reads_WeeTU.csv")
totals <- read_csv("Data/SequencedData/Plants/ProcessedData/ITS2_totals.csv")
samples <- read_csv("Data/CollectionData/fecal_sample_collection.csv")

# colorblind friendly palette
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# DATA PREP #

# get plant OTUs only
reads <- filter(reads, WTU.clade1 == 4) %>% 
  select(OTU:DataFrame)

# only PPs
samples_PP <- samples %>% filter(species == 'PP')


# Plot 1 #
# OTUs, 2017 ------------------------------------------------------------------#

# outliers are proving to be confusing

dat1 <- prep_2016_allsp_relabund_ITS2(samples, reads, totals, 1000, 2017, 0.01)
dat2 <- prep_2017_allsp_relabund_ITS2(samples, reads, totals, 1000, 2017, 0.05)
dat3 <- prep_2017_allsp_relabund_ITS2(samples, reads, totals, 1000, 2017, 0.005)
dat4 <- prep_2017_allsp_relabund_ITS2(samples, reads, totals, 1000, 2017, 0.001)
dat5 <- prep_2017_allsp_relabund_ITS2(samples, reads, totals, 2000, 2017, 0.01)
dat6 <- prep_2017_allsp_relabund_ITS2(samples, reads, totals, 2000, 2017, 0.05)
dat7 <- prep_2017_allsp_relabund_ITS2(samples, reads, totals, 2000, 2017, 0.005)
dat8 <- prep_2017_allsp_relabund_ITS2(samples, reads, totals, 2000, 2017, 0.001)
dat9 <- prep_2017_allsp_relabund_ITS2(samples, reads, totals, 5000, 2017, 0.01)
dat10 <- prep_2017_allsp_relabund_ITS2(samples, reads, totals, 5000, 2017, 0.05)
dat11 <- prep_2017_allsp_relabund_ITS2(samples, reads, totals, 5000, 2017, 0.005)
dat12 <- prep_2017_allsp_relabund_ITS2(samples, reads, totals, 5000, 2017, 0.001)

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
  ggtitle("ITS2: Spring 2017") +
  theme_bw() +
  theme(legend.position = 'bottom',
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

#ggsave("Plots/trnL_2017_allsp_totalreads_relabund.png", plot1, device = "png")

# PLOT 2 #
# OTUs, 2016 ------------------------------------------------------------------#

dat1 <- prep_2016_allsp_relabund_ITS2(samples, reads, totals, 1000, 2016, 0.025)
dat2 <- prep_2016_allsp_relabund_ITS2(samples, reads, totals, 1000, 2016, 0.01)
dat3 <- prep_2016_allsp_relabund_ITS2(samples, reads, totals, 1000, 2016, 0.005)
dat4 <- prep_2016_allsp_relabund_ITS2(samples, reads, totals, 1000, 2016, 0.001)
dat5 <- prep_2016_allsp_relabund_ITS2(samples, reads, totals, 2000, 2016, 0.01)
dat6 <- prep_2016_allsp_relabund_ITS2(samples, reads, totals, 2000, 2016, 0.025)
dat7 <- prep_2016_allsp_relabund_ITS2(samples, reads, totals, 2000, 2016, 0.005)
dat8 <- prep_2016_allsp_relabund_ITS2(samples, reads, totals, 2000, 2016, 0.001)
dat9 <- prep_2016_allsp_relabund_ITS2(samples, reads, totals, 5000, 2016, 0.01)
dat10 <- prep_2016_allsp_relabund_ITS2(samples, reads, totals, 5000, 2016, 0.025)
dat11 <- prep_2016_allsp_relabund_ITS2(samples, reads, totals, 5000, 2016, 0.005)
dat12 <- prep_2016_allsp_relabund_ITS2(samples, reads, totals, 5000, 2016, 0.001)

df <- bind_rows(dat1, dat2, dat3, dat4, dat5, dat6,
                dat7, dat8, dat9, dat10, dat11, dat12)

(plot2 <- ggplot(data = df[df$df == 'NMDS',], aes(x = MDS1, y = MDS2)) + 
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
    ggtitle("ITS2: Fall 2016") +
    theme_bw() +
    theme(legend.position = 'bottom',
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())) 

#ggsave("Plots/ITS2_2016_allsp_totalreads_relabund.png", plot2, device = "png")


# PLOT 3 #
# OTUs, 2016, PPs only ------------------------------------------------------------------#

dat1 <- prep_2016_PPonly_relabund_ITS2(samples_PP, reads, totals, 1000, 2016, 0.01)
dat2 <- prep_2016_PPonly_relabund_ITS2(samples_PP, reads, totals, 1000, 2016, 0.05)
dat3 <- prep_2016_PPonly_relabund_ITS2(samples_PP, reads, totals, 1000, 2016, 0.005)
dat4 <- prep_2016_PPonly_relabund_ITS2(samples_PP, reads, totals, 1000, 2016, 0.001)
dat5 <- prep_2016_PPonly_relabund_ITS2(samples_PP, reads, totals, 2000, 2016, 0.01)
dat6 <- prep_2016_PPonly_relabund_ITS2(samples_PP, reads, totals, 2000, 2016, 0.05)
dat7 <- prep_2016_PPonly_relabund_ITS2(samples_PP, reads, totals, 2000, 2016, 0.005)
dat8 <- prep_2016_PPonly_relabund_ITS2(samples_PP, reads, totals, 2000, 2016, 0.001)
dat9 <- prep_2016_PPonly_relabund_ITS2(samples_PP, reads, totals, 5000, 2016, 0.01)
dat10 <- prep_2016_PPonly_relabund_ITS2(samples_PP, reads, totals, 5000, 2016, 0.05)
dat11 <- prep_2016_PPonly_relabund_ITS2(samples_PP, reads, totals, 5000, 2016, 0.005)
dat12 <- prep_2016_PPonly_relabund_ITS2(samples_PP, reads, totals, 5000, 2016, 0.001)

df <- bind_rows(dat1, dat2, dat3, dat4, dat5, dat6,
                dat7, dat8, dat9, dat10, dat11, dat12)

(plot3 <- ggplot(data = df[df$df == 'NMDS',], aes(x = MDS1, y = MDS2)) + 
    geom_point(aes(color = group), size = 0.5) +
    geom_path(data = df[df$df == "df_ell",], aes(x = MDS1, y = MDS2, colour = group), 
              size = 0.5) +
    facet_wrap(min_total ~ min_rel_abund, scales = "free") +
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
    ggtitle("ITS2: Fall 2016, PP only") +
    theme_bw() +
    theme(legend.position = 'bottom',
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())) 

#ggsave("Plots/trnL_2016_allsp_totalreads_relabund.png", plot3, device = "png")

# PLOT 4 #
# OTUs, 2017, PPs only ------------------------------------------------------------------#

dat1 <- prep_2017_PPonly_relabund_ITS2(samples_PP, reads, totals, 1000, 2017, 0.01)
dat2 <- prep_2017_PPonly_relabund_ITS2(samples_PP, reads, totals, 1000, 2017, 0.05)
dat3 <- prep_2017_PPonly_relabund_ITS2(samples_PP, reads, totals, 1000, 2017, 0.005)
dat4 <- prep_2017_PPonly_relabund_ITS2(samples_PP, reads, totals, 1000, 2017, 0.001)
dat5 <- prep_2017_PPonly_relabund_ITS2(samples_PP, reads, totals, 2000, 2017, 0.01)
dat6 <- prep_2017_PPonly_relabund_ITS2(samples_PP, reads, totals, 2000, 2017, 0.05)
dat7 <- prep_2017_PPonly_relabund_ITS2(samples_PP, reads, totals, 2000, 2017, 0.005)
dat8 <- prep_2017_PPonly_relabund_ITS2(samples_PP, reads, totals, 2000, 2017, 0.001)
dat9 <- prep_2017_PPonly_relabund_ITS2(samples_PP, reads, totals, 5000, 2017, 0.01)
dat10 <- prep_2017_PPonly_relabund_ITS2(samples_PP, reads, totals, 5000, 2017, 0.05)
dat11 <- prep_2017_PPonly_relabund_ITS2(samples_PP, reads, totals, 5000, 2017, 0.005)
dat12 <- prep_2017_PPonly_relabund_ITS2(samples_PP, reads, totals, 5000, 2017, 0.001)

df <- bind_rows(dat1, dat2, dat3, dat4, dat5, dat6,
                dat7, dat8, dat9, dat10, dat11, dat12)

(plot4 <- ggplot(data = df[df$df == 'NMDS',], aes(x = MDS1, y = MDS2)) + 
    geom_point(aes(color = group), size = 0.5) +
    geom_path(data = df[df$df == "df_ell",], aes(x = MDS1, y = MDS2, colour = group), 
              size = 0.5) +
    facet_wrap(min_total ~ min_rel_abund, scales = "free") +
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
    ggtitle("ITS2: Spring 2017, PP only") +
    theme_bw() +
    theme(legend.position = 'bottom',
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())) 

#ggsave("Plots/trnL_2017_PPonly_totalreads_relabund.png", plot4, device = "png")

# WORKING AREA ================================================================#

# for finding outliers
# 
# data <- filter_reads_data_ITS2(samples,
#                                reads,
#                                totals,
#                                reads_min = 5000,
#                                yr = 2016,
#                                rel_reads_min = 0.05) %>%
#   data_prep_multivariate()
# data[[1]] <- binarize(data[[1]])
# 
# # remove outliers
# # "S008810", "S010014"
# data[[1]] <- 
#   data[[1]][!(row.names(data[[1]]) %in% c("S008824")),]
# data[[2]] <- 
#   data[[2]][!data[[2]] %in% c("S008824")]
# data[[3]] <- 
#   data[[3]][!(data[[3]]$vial_barcode) %in% c("S008824"),]
# 
# dist_trnL <- metaMDS(data[[1]], distance = "bray", trymax = 250, k = 3)
# dist_trnL$points
# 
# # scree plot
# goeveg::dimcheckMDS(data[[1]], distance = "bray", k = 6, trymax = 50)
