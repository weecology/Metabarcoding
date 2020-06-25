# NMDS trnL Markdown Plotting
# EKB
# May 2020

# LIBRARIES and DATA #

library(tidyverse)
library(vegan)
source('Scripts/functions.R')

reads <- read_csv("Data/SequencedData/Plants/ProcessedData/trnL_reads.csv")
totals <- read_csv("Data/SequencedData/Plants/ProcessedData/trnL_totals.csv")
samples <- read_csv("Data/CollectionData/fecal_sample_collection.csv")

# colorblind friendly palette
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# DATA PREP #

# only PPs
samples_PP <- samples %>% filter(species == 'PP')

# OTUS #========================================================================

# Plot 1 #
# OTUs, 2017 ------------------------------------------------------------------#

dat1 <- prep_2017_allsp_relabund(samples, reads, totals, 1000, 2017, 0.01)
dat2 <- prep_2017_allsp_relabund(samples, reads, totals, 1000, 2017, 0.05)
dat3 <- prep_2017_allsp_relabund(samples, reads, totals, 1000, 2017, 0.005)
dat4 <- prep_2017_allsp_relabund(samples, reads, totals, 1000, 2017, 0.001)
dat5 <- prep_2017_allsp_relabund(samples, reads, totals, 2000, 2017, 0.01)
dat6 <- prep_2017_allsp_relabund(samples, reads, totals, 2000, 2017, 0.05)
dat7 <- prep_2017_allsp_relabund(samples, reads, totals, 2000, 2017, 0.005)
dat8 <- prep_2017_allsp_relabund(samples, reads, totals, 2000, 2017, 0.001)
dat9 <- prep_2017_allsp_relabund(samples, reads, totals, 5000, 2017, 0.01)
dat10 <- prep_2017_allsp_relabund(samples, reads, totals, 5000, 2017, 0.05)
dat11 <- prep_2017_allsp_relabund(samples, reads, totals, 5000, 2017, 0.005)
dat12 <- prep_2017_allsp_relabund(samples, reads, totals, 5000, 2017, 0.001)

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
  ggtitle("trnL: Spring 2017") +
  theme_bw() +
  theme(legend.position = 'bottom',
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

#ggsave("Plots/trnL_2017_allsp_totalreads_relabund.png", plot1, device = "png")

# PLOT 2 #
# OTUs, 2016 ------------------------------------------------------------------#

dat1 <- prep_2016_allsp_relabund(samples, reads, totals, 1000, 2016, 0.01)
dat2 <- prep_2016_allsp_relabund(samples, reads, totals, 1000, 2016, 0.05)
dat3 <- prep_2016_allsp_relabund(samples, reads, totals, 1000, 2016, 0.005)
dat4 <- prep_2016_allsp_relabund(samples, reads, totals, 1000, 2016, 0.001)
dat5 <- prep_2016_allsp_relabund(samples, reads, totals, 2000, 2016, 0.01)
dat6 <- prep_2016_allsp_relabund(samples, reads, totals, 2000, 2016, 0.05)
dat7 <- prep_2016_allsp_relabund(samples, reads, totals, 2000, 2016, 0.005)
dat8 <- prep_2016_allsp_relabund(samples, reads, totals, 2000, 2016, 0.001)
dat9 <- prep_2016_allsp_relabund(samples, reads, totals, 5000, 2016, 0.01)
dat10 <- prep_2016_allsp_relabund(samples, reads, totals, 5000, 2016, 0.05)
dat11 <- prep_2016_allsp_relabund(samples, reads, totals, 5000, 2016, 0.005)
dat12 <- prep_2016_allsp_relabund(samples, reads, totals, 5000, 2016, 0.001)

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
  ggtitle("trnL: Fall 2016") +
  theme_bw() +
  theme(legend.position = 'bottom',
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())) 

#ggsave("Plots/trnL_2016_allsp_totalreads_relabund.png", plot2, device = "png")


# PLOT 3 #
# OTUs, 2016, PPs only ------------------------------------------------------------------#

dat1 <- prep_2016_PPonly_relabund(samples_PP, reads, totals, 1000, 2016, 0.01)
dat2 <- prep_2016_PPonly_relabund(samples_PP, reads, totals, 1000, 2016, 0.05)
dat3 <- prep_2016_PPonly_relabund(samples_PP, reads, totals, 1000, 2016, 0.005)
dat4 <- prep_2016_PPonly_relabund(samples_PP, reads, totals, 1000, 2016, 0.001)
dat5 <- prep_2016_PPonly_relabund(samples_PP, reads, totals, 2000, 2016, 0.01)
dat6 <- prep_2016_PPonly_relabund(samples_PP, reads, totals, 2000, 2016, 0.05)
dat7 <- prep_2016_PPonly_relabund(samples_PP, reads, totals, 2000, 2016, 0.005)
dat8 <- prep_2016_PPonly_relabund(samples_PP, reads, totals, 2000, 2016, 0.001)
dat9 <- prep_2016_PPonly_relabund(samples_PP, reads, totals, 5000, 2016, 0.01)
dat10 <- prep_2016_PPonly_relabund(samples_PP, reads, totals, 5000, 2016, 0.05)
dat11 <- prep_2016_PPonly_relabund(samples_PP, reads, totals, 5000, 2016, 0.005)
dat12 <- prep_2016_PPonly_relabund(samples_PP, reads, totals, 5000, 2016, 0.001)

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
    ggtitle("trnL: Fall 2016, PP only") +
    theme_bw() +
    theme(legend.position = 'bottom',
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())) 

#ggsave("Plots/trnL_2016_PPonly_totalreads_relabund.png", plot3, device = "png")

# PLOT 4 #
# OTUs, 2017, PPs only ------------------------------------------------------------------#

dat1 <- prep_2017_PPonly_relabund(samples_PP, reads, totals, 1000, 2017, 0.01)
dat2 <- prep_2017_PPonly_relabund(samples_PP, reads, totals, 1000, 2017, 0.05)
dat3 <- prep_2017_PPonly_relabund(samples_PP, reads, totals, 1000, 2017, 0.005)
dat4 <- prep_2017_PPonly_relabund(samples_PP, reads, totals, 1000, 2017, 0.001)
dat5 <- prep_2017_PPonly_relabund(samples_PP, reads, totals, 2000, 2017, 0.01)
dat6 <- prep_2017_PPonly_relabund(samples_PP, reads, totals, 2000, 2017, 0.05)
dat7 <- prep_2017_PPonly_relabund(samples_PP, reads, totals, 2000, 2017, 0.005)
dat8 <- prep_2017_PPonly_relabund(samples_PP, reads, totals, 2000, 2017, 0.001)
dat9 <- prep_2017_PPonly_relabund(samples_PP, reads, totals, 5000, 2017, 0.01)
dat10 <- prep_2017_PPonly_relabund(samples_PP, reads, totals, 5000, 2017, 0.05)
dat11 <- prep_2017_PPonly_relabund(samples_PP, reads, totals, 5000, 2017, 0.005)
dat12 <- prep_2017_PPonly_relabund(samples_PP, reads, totals, 5000, 2017, 0.001)

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
    ggtitle("trnL: Spring 2017, PP only") +
    theme_bw() +
    theme(legend.position = 'bottom',
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())) 

#ggsave("Plots/trnL_2017_PPonly_totalreads_relabund.png", plot4, device = "png")
