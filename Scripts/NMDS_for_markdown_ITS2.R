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


# OTUS #========================================================================

# PLOT 1 #
# OTUs, 2016 ------------------------------------------------------------------#

dat1 <- prep_454_allsp_relabund_ITS2(samples, reads, totals, 1000, 454, 0.01)
dat2 <- prep_454_allsp_relabund_ITS2(samples, reads, totals, 1000, 454, 0.05)
dat3 <- prep_454_allsp_relabund_ITS2(samples, reads, totals, 1000, 454, 0.005)
dat4 <- prep_454_allsp_relabund_ITS2(samples, reads, totals, 1000, 454, 0.001)
dat5 <- prep_454_allsp_relabund_ITS2(samples, reads, totals, 2000, 454, 0.01)
dat6 <- prep_454_allsp_relabund_ITS2(samples, reads, totals, 2000, 454, 0.05)
dat7 <- prep_454_allsp_relabund_ITS2(samples, reads, totals, 2000, 454, 0.005)
dat8 <- prep_454_allsp_relabund_ITS2(samples, reads, totals, 2000, 454, 0.001)
dat9 <- prep_454_allsp_relabund_ITS2(samples, reads, totals, 5000, 454, 0.01)
dat10 <- prep_454_allsp_relabund_ITS2(samples, reads, totals, 5000, 454, 0.05)
dat11 <- prep_454_allsp_relabund_ITS2(samples, reads, totals, 5000, 454, 0.005)
dat12 <- prep_454_allsp_relabund_ITS2(samples, reads, totals, 5000, 454, 0.001)

df <- bind_rows(dat1, dat2, dat3, dat4, dat5, dat6,
                dat7, dat8, dat9, dat10, dat11, dat12)

(plot1 <- ggplot(data = df[df$df == 'NMDS',], aes(x = MDS1, y = MDS2)) + 
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

#ggsave("Plots/ITS2_454_allsp_totalreads_relabund.png", plot1, device = "png")

# PLOT 2 #
# OTUs, 454, PPs only ------------------------------------------------------------------#

dat1 <- prep_454_PPonly_relabund_ITS2(samples_PP, reads, totals, 1000, 454, 0.01)
dat3 <- prep_454_PPonly_relabund_ITS2(samples_PP, reads, totals, 1000, 454, 0.005)
dat4 <- prep_454_PPonly_relabund_ITS2(samples_PP, reads, totals, 1000, 454, 0.001)
dat5 <- prep_454_PPonly_relabund_ITS2(samples_PP, reads, totals, 2000, 454, 0.01)
dat7 <- prep_454_PPonly_relabund_ITS2(samples_PP, reads, totals, 2000, 454, 0.005)
dat8 <- prep_454_PPonly_relabund_ITS2(samples_PP, reads, totals, 2000, 454, 0.001)
dat9 <- prep_454_PPonly_relabund_ITS2(samples_PP, reads, totals, 5000, 454, 0.01)
dat11 <- prep_454_PPonly_relabund_ITS2(samples_PP, reads, totals, 5000, 454, 0.005)
dat12 <- prep_454_PPonly_relabund_ITS2(samples_PP, reads, totals, 5000, 454, 0.001)

df <- bind_rows(dat1, dat3, dat4, dat5,
                dat7, dat8, dat9, dat11, dat12)

(plot2 <- ggplot(data = df[df$df == 'NMDS',], aes(x = MDS1, y = MDS2)) + 
    geom_point(aes(color = group), size = 0.5) +
    geom_path(data = df[df$df == "df_ell",], aes(x = MDS1, y = MDS2, colour = group), 
              size = 0.5) +
    facet_wrap(min_total ~ min_rel_abund, scales = "free") +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "PP: Control",], 
              aes(x = MDS1, y = MDS2,  
                  label = .data$group[1], 
                  color = .data$group[1]),
              size = 1) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "PP: KR_exclosure",],
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

#ggsave("Plots/ITS2_2016_PPonly_totalreads_relabund.png", plot2, device = "png")

# Plot 3 #
# OTUs, 460 ------------------------------------------------------------------#

dat1 <- prep_460_allsp_relabund_ITS2(samples, reads, totals, 1000, 460, 0.01)
dat3 <- prep_460_allsp_relabund_ITS2(samples, reads, totals, 1000, 460, 0.005)
dat4 <- prep_460_allsp_relabund_ITS2(samples, reads, totals, 1000, 460, 0.001)
dat5 <- prep_460_allsp_relabund_ITS2(samples, reads, totals, 2000, 460, 0.01)
dat7 <- prep_460_allsp_relabund_ITS2(samples, reads, totals, 2000, 460, 0.005)
dat8 <- prep_460_allsp_relabund_ITS2(samples, reads, totals, 2000, 460, 0.001)
dat9 <- prep_460_allsp_relabund_ITS2(samples, reads, totals, 5000, 460, 0.01)
dat11 <- prep_460_allsp_relabund_ITS2(samples, reads, totals, 5000, 460, 0.005)
dat12 <- prep_460_allsp_relabund_ITS2(samples, reads, totals, 5000, 460, 0.001)

df<- bind_rows(dat1, dat3, dat4, dat5,
               dat7, dat8, dat9, dat11, dat12)

(plot3 <- ggplot(data = df[df$df == 'NMDS',], aes(x = MDS1, y = MDS2)) + 
    geom_point(aes(color = group), size = 0.5) +
    geom_path(data = df[df$df == "df_ell",], aes(x = MDS1, y = MDS2, colour = group), 
              size = 0.5) +
    facet_wrap(min_total ~ min_rel_abund, scales = "free", ncol = 3) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "K-Rat",], 
              aes(x = MDS1, y = MDS2,  
                  label = .data$group[1], 
                  color = .data$group[1]),
              size = 1) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "PP: Control",], 
              aes(x = MDS1, y = MDS2,
                  label = .data$group[1], 
                  color = .data$group[1]),
              size = 1) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "PP: KR_Exclosure",],
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
          panel.grid.minor = element_blank())) 

#ggsave("Plots/ITS2_460_allsp_totalreads_relabund.png", plot3, device = "png")

# PLOT 4 #
# OTUs, 460, PPs only ------------------------------------------------------------------#

dat1 <- prep_460_PPonly_relabund_ITS2(samples_PP, reads, totals, 1000, 460, 0.01)
dat3 <- prep_460_PPonly_relabund_ITS2(samples_PP, reads, totals, 1000, 460, 0.005)
dat4 <- prep_460_PPonly_relabund_ITS2(samples_PP, reads, totals, 1000, 460, 0.001)
dat5 <- prep_460_PPonly_relabund_ITS2(samples_PP, reads, totals, 2000, 460, 0.01)
dat7 <- prep_460_PPonly_relabund_ITS2(samples_PP, reads, totals, 2000, 460, 0.005)
dat8 <- prep_460_PPonly_relabund_ITS2(samples_PP, reads, totals, 2000, 460, 0.001)
dat9 <- prep_460_PPonly_relabund_ITS2(samples_PP, reads, totals, 5000, 460, 0.01)
dat11 <- prep_460_PPonly_relabund_ITS2(samples_PP, reads, totals, 5000, 460, 0.005)
dat12 <- prep_460_PPonly_relabund_ITS2(samples_PP, reads, totals, 5000, 460, 0.001)

df <- bind_rows(dat1, dat3, dat4, dat5,
                dat7, dat8, dat9, dat11, dat12)

(plot4 <- ggplot(data = df[df$df == 'NMDS',], aes(x = MDS1, y = MDS2)) + 
    geom_point(aes(color = group), size = 0.5) +
    geom_path(data = df[df$df == "df_ell",], aes(x = MDS1, y = MDS2, colour = group), 
              size = 0.5) +
    facet_wrap(min_total ~ min_rel_abund, scales = "free") +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "PP: Control",], 
              aes(x = MDS1, y = MDS2,  
                  label = .data$group[1], 
                  color = .data$group[1]),
              size = 1) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "PP: KR_Exclosure",],
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

#ggsave("Plots/ITS2_460_PPonly_totalreads_relabund.png", plot4, device = "png")

# Plot 5 #
# OTUs, 466 ------------------------------------------------------------------#

dat1 <- prep_466_allsp_relabund_ITS2(samples, reads, totals, 1000, 466, 0.01)
dat3 <- prep_466_allsp_relabund_ITS2(samples, reads, totals, 1000, 466, 0.005)
dat4 <- prep_466_allsp_relabund_ITS2(samples, reads, totals, 1000, 466, 0.001)
dat5 <- prep_466_allsp_relabund_ITS2(samples, reads, totals, 2000, 466, 0.01)
dat7 <- prep_466_allsp_relabund_ITS2(samples, reads, totals, 2000, 466, 0.005)
dat8 <- prep_466_allsp_relabund_ITS2(samples, reads, totals, 2000, 466, 0.001)
dat9 <- prep_466_allsp_relabund_ITS2(samples, reads, totals, 5000, 466, 0.01)
dat11 <- prep_466_allsp_relabund_ITS2(samples, reads, totals, 5000, 466, 0.005)
dat12 <- prep_466_allsp_relabund_ITS2(samples, reads, totals, 5000, 466, 0.001)

df<- bind_rows(dat1, dat3, dat4, dat5,
               dat7, dat8, dat9, dat11, dat12)

(plot5 <- ggplot(data = df[df$df == 'NMDS',], aes(x = MDS1, y = MDS2)) + 
    geom_point(aes(color = group), size = 0.5) +
    geom_path(data = df[df$df == "df_ell",], aes(x = MDS1, y = MDS2, colour = group), 
              size = 0.5) +
    facet_wrap(min_total ~ min_rel_abund, scales = "free", ncol = 3) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "K-Rat",], 
              aes(x = MDS1, y = MDS2,  
                  label = .data$group[1], 
                  color = .data$group[1]),
              size = 1) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "PP: Control",], 
              aes(x = MDS1, y = MDS2,
                  label = .data$group[1], 
                  color = .data$group[1]),
              size = 1) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "PP: KR_Exclosure",],
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
    ggtitle("ITS2: Fall 2017") +
    theme_bw() +
    theme(legend.position = 'bottom',
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())) 

#ggsave("Plots/ITS2_466_allsp_totalreads_relabund.png", plot5, device = "png")

# PLOT 6#
# OTUs, 466, PPs only ------------------------------------------------------------------#

dat1 <- prep_466_PPonly_relabund_ITS2(samples_PP, reads, totals, 1000, 466, 0.01)
dat3 <- prep_466_PPonly_relabund_ITS2(samples_PP, reads, totals, 1000, 466, 0.005)
dat4 <- prep_466_PPonly_relabund_ITS2(samples_PP, reads, totals, 1000, 466, 0.001)
dat5 <- prep_466_PPonly_relabund_ITS2(samples_PP, reads, totals, 2000, 466, 0.01)
dat7 <- prep_466_PPonly_relabund_ITS2(samples_PP, reads, totals, 2000, 466, 0.005)
dat8 <- prep_466_PPonly_relabund_ITS2(samples_PP, reads, totals, 2000, 466, 0.001)
dat9 <- prep_466_PPonly_relabund_ITS2(samples_PP, reads, totals, 5000, 466, 0.01)
dat11 <- prep_466_PPonly_relabund_ITS2(samples_PP, reads, totals, 5000, 466, 0.005)
dat12 <- prep_466_PPonly_relabund_ITS2(samples_PP, reads, totals, 5000, 466, 0.001)

df <- bind_rows(dat1, dat3, dat4, dat5,
                dat7, dat8, dat9, dat11, dat12)

(plot6 <- ggplot(data = df[df$df == 'NMDS',], aes(x = MDS1, y = MDS2)) + 
    geom_point(aes(color = group), size = 0.5) +
    geom_path(data = df[df$df == "df_ell",], aes(x = MDS1, y = MDS2, colour = group), 
              size = 0.5) +
    facet_wrap(min_total ~ min_rel_abund, scales = "free") +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "PP: Control",], 
              aes(x = MDS1, y = MDS2,  
                  label = .data$group[1], 
                  color = .data$group[1]),
              size = 1) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "PP: KR_Exclosure",],
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
    ggtitle("ITS2: Fall 2017, PP only") +
    theme_bw() +
    theme(legend.position = 'bottom',
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())) 

#ggsave("Plots/ITS2_466_PPonly_totalreads_relabund.png", plot6, device = "png")


# WORKING AREA ================================================================#

# for finding outliers

data <- filter_reads_data_ITS2(samples,
                               reads,
                               totals,
                               reads_min = 2000,
                               period_code = 466,
                               rel_reads_min = 0.01) %>%
  data_prep_multivariate()
data[[1]] <- binarize(data[[1]])

# remove outliers
data[[1]] <- 
  data[[1]][!(row.names(data[[1]]) %in% c("S013043")),]
data[[2]] <- 
  data[[2]][!data[[2]] %in% c("S013043")]
data[[3]] <- 
  data[[3]][!(data[[3]]$vial_barcode) %in% c("S013043"),]

dist_trnL <- metaMDS(data[[1]], distance = "bray", trymax = 250, k = 3)
dist_matrix <- metaMDSredist(dist_trnL)
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

# run perMANOVA
group = as.matrix(groups$group)
perMANOVA_output <- adonis(data[[1]] ~ group, permutations = 10000)
pairwise_perMANOVA <- adonis.pair(dist.mat = dist_matrix, Factor = as.factor(groups$group))
pairwise_perMANOVA
