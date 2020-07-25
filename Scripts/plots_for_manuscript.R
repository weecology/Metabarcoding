# Make Figures for Dissertation
# EKB
# July/August 2020

library(tidyverse)
library(patchwork)
library(vegan)
source('Scripts/functions.R')

# colorblind friendly palette
cbPalette <- c("#009E73", "#E69F00", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# TRNL OTUS #

reads <- read_csv("Data/SequencedData/Plants/ProcessedData/trnL_reads.csv")
totals <- read_csv("Data/SequencedData/Plants/ProcessedData/trnL_totals.csv")
samples <- read_csv("Data/CollectionData/fecal_sample_collection.csv")

# only PPs
samples_PP <- samples %>% filter(species == 'PP')

# plot 1
df <- prep_2016_allsp_relabund(samples, reads, totals, 2000, 2016, 0.01)
(plot1 <- ggplot(data = df[df$df == 'NMDS',], aes(x = MDS1, y = MDS2)) + 
    geom_point(aes(color = group), size = 0.5, alpha = 0.5) +
    geom_path(data = df[df$df == "df_ell",], aes(x = MDS1, y = MDS2, colour = group), 
              size = 0.5) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "K-Rat",], 
              aes(x = MDS1, y = MDS2,  
                  label = .data$group[1], 
                  color = .data$group[1]),
              size = 2) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "PP: Control",], 
              aes(x = MDS1, y = MDS2,  
                  label = .data$group[1], 
                  color = .data$group[1]),
              size = 2) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "PP: KR_Exclosure",],
              aes(x = MDS1, y = MDS2,  
                  label = .data$group[1], 
                  color = .data$group[1]),
              size = 2) +
    geom_text(data = df %>% select(F.model, pval, min_total, min_rel_abund) %>% distinct(), 
              aes(x = Inf, y = Inf,
                  label = paste("F.model = ", round(.data$F.model, 2),
                                "\n p = ", round(.data$pval, 4))),
              hjust = 1.1, vjust= 1.2, size = 2) +
    scale_color_manual(values = cbPalette) +
    labs(tag = "A") +
    theme_bw() +
    theme(legend.title = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank())) 

legend <- cowplot::get_legend(plot1)

plot1 <- plot1 + theme(legend.position = "none")

# plot 2
df <- prep_2017_allsp_relabund(samples, reads, totals, 2000, 2017, 0.01)
(plot2 <- ggplot(data = df[df$df == 'NMDS',], aes(x = MDS1, y = MDS2)) + 
  geom_point(aes(color = group), size = 0.5, alpha = 0.5) +
  geom_path(data = df[df$df == "df_ell",], aes(x = MDS1, y = MDS2, colour = group), 
            size = 0.5) +
  geom_text(data = df[df$df == "NMDS.mean" & df$group == "K-Rat",], 
            aes(x = MDS1, y = MDS2,  
                label = .data$group[1], 
                color = .data$group[1]),
            size = 2) +
  geom_text(data = df[df$df == "NMDS.mean" & df$group == "PP: Control",], 
            aes(x = MDS1, y = MDS2,  
                label = .data$group[1], 
                color = .data$group[1]),
            size = 2) +
  geom_text(data = df[df$df == "NMDS.mean" & df$group == "PP: KR_Exclosure",],
            aes(x = MDS1, y = MDS2,  
                label = .data$group[1], 
                color = .data$group[1]),
            size = 2) +
  geom_text(data = df %>% select(F.model, pval, min_total, min_rel_abund) %>% distinct(), 
            aes(x = Inf, y = Inf,
                label = paste("F.model = ", round(.data$F.model, 2),
                              "\n p = ", round(.data$pval, 4))),
            hjust = 1.1, vjust= 1.2, size = 2) +
  scale_color_manual(values = c("#009E73", "#E69F00", "#56B4E9")) +
  labs(tag = "B") +  
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())) 

# plot 3
df <- prep_2016_PPonly_relabund(samples_PP, reads, totals, 2000, 2016, 0.01)
(plot3 <- ggplot(data = df[df$df == 'NMDS',], aes(x = MDS1, y = MDS2)) + 
    geom_point(aes(color = group), size = 0.5, alpha = 0.5) +
    geom_path(data = df[df$df == "df_ell",], aes(x = MDS1, y = MDS2, colour = group), 
              size = 0.5) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "PP: Control",], 
              aes(x = MDS1, y = MDS2,  
                  label = .data$group[1], 
                  color = .data$group[1]),
              size = 2) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "PP: KR_Exclosure",],
              aes(x = MDS1, y = MDS2,  
                  label = .data$group[1], 
                  color = .data$group[1]),
              size = 2) +
    geom_text(data = df %>% select(F.model, pval, min_total, min_rel_abund) %>% distinct(), 
              aes(x = Inf, y = Inf,
                  label = paste("F.model = ", round(.data$F.model, 2),
                                "\n p = ", round(.data$pval, 4))),
              hjust = 1.1, vjust= 1.2, size = 2) +
    scale_color_manual(values = c("#E69F00", "#56B4E9")) +
    labs(tag = "C") +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())) 

# plot 4
df <- prep_2017_PPonly_relabund(samples_PP, reads, totals, 2000, 2017, 0.01)
(plot4 <- ggplot(data = df[df$df == 'NMDS',], aes(x = MDS1, y = MDS2)) + 
    geom_point(aes(color = group), size = 0.5, alpha = 0.5) +
    geom_path(data = df[df$df == "df_ell",], aes(x = MDS1, y = MDS2, colour = group), 
              size = 0.5) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "PP: Control",], 
              aes(x = MDS1, y = MDS2,  
                  label = .data$group[1], 
                  color = .data$group[1]),
              size = 2) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "PP: KR_Exclosure",],
              aes(x = MDS1, y = MDS2,  
                  label = .data$group[1], 
                  color = .data$group[1]),
              size = 2) +
    geom_text(data = df %>% select(F.model, pval, min_total, min_rel_abund) %>% distinct(), 
              aes(x = Inf, y = Inf,
                  label = paste("F.model = ", round(.data$F.model, 2),
                                "\n p = ", round(.data$pval, 4))),
              hjust = 1.1, vjust= 1.2, size = 2) +
    scale_color_manual(values = c("#E69F00", "#56B4E9")) +
    labs(tag = "D") +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank())) 

# multipanel plot

col1 <- ggplot() + annotate(geom = 'text', x = 1, y = 1, 
                            label = "Fall 2016", size = 5, fontface = "bold") + 
  theme_void() 
col2 <- ggplot() + annotate(geom = 'text', x = 1, y = 1, 
                            label = "Spring 2017", size = 5, fontface = "bold") + 
  theme_void() 

layoutplot <- "
cccddd##
eeefff##
eeefffjj
ggghhhjj
ggghhh##
"
plot1 <- plot1 + theme(plot.margin = unit(c(0,20,0,0), "pt"))

plotlist <- list(c = col1, d = col2, e = plot1, f = plot2, 
                 g = plot3, h = plot4, j = legend)
wrap_plots(plotlist, design = layoutplot) 

#ggsave("Plots/trnL_OTUs.png")

# ITS2 OTUS #-------------------------------------------------------------------

reads <- read_csv("Data/SequencedData/Plants/ProcessedData/ITS2_reads_WeeTU.csv")
totals <- read_csv("Data/SequencedData/Plants/ProcessedData/ITS2_totals.csv")
samples <- read_csv("Data/CollectionData/fecal_sample_collection.csv")

# get plant OTUs only
reads <- filter(reads, WTU.clade1 == 4) %>% 
  select(OTU:DataFrame)

# only PPs
samples_PP <- samples %>% filter(species == 'PP')

# plot 1
df <- prep_2016_allsp_relabund_ITS2(samples, reads, totals, 2000, 2016, 0.01)
(plot1 <- ggplot(data = df[df$df == 'NMDS',], aes(x = MDS1, y = MDS2)) + 
    geom_point(aes(color = group), size = 0.5, alpha = 0.5) +
    geom_path(data = df[df$df == "df_ell",], aes(x = MDS1, y = MDS2, colour = group), 
              size = 0.5) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "K-Rat",], 
              aes(x = MDS1, y = MDS2,  
                  label = .data$group[1], 
                  color = .data$group[1]),
              size = 2) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "PP: Control",], 
              aes(x = MDS1, y = MDS2,  
                  label = .data$group[1], 
                  color = .data$group[1]),
              size = 2) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "PP: KR_Exclosure",],
              aes(x = MDS1, y = MDS2,  
                  label = .data$group[1], 
                  color = .data$group[1]),
              size = 2) +
    geom_text(data = df %>% select(F.model, pval, min_total, min_rel_abund) %>% distinct(), 
              aes(x = Inf, y = Inf,
                  label = paste("F.model = ", round(.data$F.model, 2),
                                "\n p = ", round(.data$pval, 4))),
              hjust = 1.1, vjust= 1.2, size = 2) +
    scale_color_manual(values = cbPalette) +
    labs(tag = "A") +
    theme_bw() +
    theme(legend.title = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank())) 

legend <- cowplot::get_legend(plot1)

plot1 <- plot1 + theme(legend.position = "none")

# plot 2
df <- prep_2017_allsp_relabund_ITS2(samples, reads, totals, 2000, 2017, 0.01)
(plot2 <- ggplot(data = df[df$df == 'NMDS',], aes(x = MDS1, y = MDS2)) + 
    geom_point(aes(color = group), size = 0.5, alpha = 0.5) +
    geom_path(data = df[df$df == "df_ell",], aes(x = MDS1, y = MDS2, colour = group), 
              size = 0.5) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "K-Rat",], 
              aes(x = MDS1, y = MDS2,  
                  label = .data$group[1], 
                  color = .data$group[1]),
              size = 2) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "PP: Control",], 
              aes(x = MDS1, y = MDS2,  
                  label = .data$group[1], 
                  color = .data$group[1]),
              size = 2) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "PP: KR_Exclosure",],
              aes(x = MDS1, y = MDS2,  
                  label = .data$group[1], 
                  color = .data$group[1]),
              size = 2) +
    geom_text(data = df %>% select(F.model, pval, min_total, min_rel_abund) %>% distinct(), 
              aes(x = Inf, y = Inf,
                  label = paste("F.model = ", round(.data$F.model, 2),
                                "\n p = ", round(.data$pval, 4))),
              hjust = 1.1, vjust= 1.2, size = 2) +
    scale_color_manual(values = c("#009E73", "#E69F00", "#56B4E9")) +
    labs(tag = "B") +  
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank())) 

# plot 3
df <- prep_2016_PPonly_relabund_ITS2(samples_PP, reads, totals, 2000, 2016, 0.01)
(plot3 <- ggplot(data = df[df$df == 'NMDS',], aes(x = MDS1, y = MDS2)) + 
    geom_point(aes(color = group), size = 0.5, alpha = 0.5) +
    geom_path(data = df[df$df == "df_ell",], aes(x = MDS1, y = MDS2, colour = group), 
              size = 0.5) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "PP: Control",], 
              aes(x = MDS1, y = MDS2,  
                  label = .data$group[1], 
                  color = .data$group[1]),
              size = 2) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "PP: KR_Exclosure",],
              aes(x = MDS1, y = MDS2,  
                  label = .data$group[1], 
                  color = .data$group[1]),
              size = 2) +
    geom_text(data = df %>% select(F.model, pval, min_total, min_rel_abund) %>% distinct(), 
              aes(x = Inf, y = Inf,
                  label = paste("F.model = ", round(.data$F.model, 2),
                                "\n p = ", round(.data$pval, 4))),
              hjust = 1.1, vjust= 1.2, size = 2) +
    scale_color_manual(values = c("#E69F00", "#56B4E9")) +
    labs(tag = "C") +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())) 

# plot 4
df <- prep_2017_PPonly_relabund_ITS2(samples_PP, reads, totals, 2000, 2017, 0.01)
(plot4 <- ggplot(data = df[df$df == 'NMDS',], aes(x = MDS1, y = MDS2)) + 
    geom_point(aes(color = group), size = 0.5, alpha = 0.5) +
    geom_path(data = df[df$df == "df_ell",], aes(x = MDS1, y = MDS2, colour = group), 
              size = 0.5) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "PP: Control",], 
              aes(x = MDS1, y = MDS2,  
                  label = .data$group[1], 
                  color = .data$group[1]),
              size = 2) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "PP: KR_Exclosure",],
              aes(x = MDS1, y = MDS2,  
                  label = .data$group[1], 
                  color = .data$group[1]),
              size = 2) +
    geom_text(data = df %>% select(F.model, pval, min_total, min_rel_abund) %>% distinct(), 
              aes(x = Inf, y = Inf,
                  label = paste("F.model = ", round(.data$F.model, 2),
                                "\n p = ", round(.data$pval, 4))),
              hjust = 1.1, vjust= 1.2, size = 2) +
    scale_color_manual(values = c("#E69F00", "#56B4E9")) +
    labs(tag = "D") +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank())) 

# multipanel plot

col1 <- ggplot() + annotate(geom = 'text', x = 1, y = 1, 
                            label = "Fall 2016", size = 5, fontface = "bold") + 
  theme_void() 
col2 <- ggplot() + annotate(geom = 'text', x = 1, y = 1, 
                            label = "Spring 2017", size = 5, fontface = "bold") + 
  theme_void() 

layoutplot <- "
cccddd##
eeefff##
eeefffjj
ggghhhjj
ggghhh##
"
plot1 <- plot1 + theme(plot.margin = unit(c(0,20,0,0), "pt"))

plotlist <- list(c = col1, d = col2, e = plot1, f = plot2, 
                 g = plot3, h = plot4, j = legend)
wrap_plots(plotlist, design = layoutplot) 

#ggsave("Plots/ITS2_OTUs.png")
