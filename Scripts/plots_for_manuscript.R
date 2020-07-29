# Make Figures for Dissertation
# EKB
# July/August 2020

library(tidyverse)
library(patchwork)
library(vegan)
library(EcolUtils)
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
list1 <- prep_454_allsp_relabund(samples, reads, totals, 2000, 454, 0.01)
df <- list1[[1]]
(plot1 <- ggplot(data = df[df$df == 'NMDS',], aes(x = MDS1, y = MDS2)) + 
    geom_point(aes(color = group), size = 0.5, alpha = 0.5) +
    geom_path(data = df[df$df == "df_ell",], aes(x = MDS1, y = MDS2, colour = group), 
              size = 0.5) +
    # geom_text(data = df[df$df == "NMDS.mean" & df$group == "K-Rat",], 
    #           aes(x = MDS1, y = MDS2,  
    #               label = .data$group[1], 
    #               color = .data$group[1]),
    #           size = 2) +
    # geom_text(data = df[df$df == "NMDS.mean" & df$group == "CP: Control",], 
    #           aes(x = MDS1, y = MDS2,  
    #               label = .data$group[1], 
    #               color = .data$group[1]),
    #           size = 2) +
    # geom_text(data = df[df$df == "NMDS.mean" & df$group == "CP: KR Exclosure",],
    #           aes(x = MDS1, y = MDS2,  
    #               label = .data$group[1], 
    #               color = .data$group[1]),
    #           size = 2) +
    geom_text(data = df %>% select(F.model, pval, min_total, min_rel_abund) %>% distinct(), 
              aes(x = Inf, y = Inf,
                  label = paste("F.model = ", round(.data$F.model, 2),
                                "\n p = ", round(.data$pval, 4))),
              hjust = 1.1, vjust= 1.2, size = 2) +
    scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
    labs(tag = "A") +
    theme_bw() +
    theme(legend.title = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank())) 

legend <- cowplot::get_legend(plot1)
plot1 <- plot1 + theme(legend.position = "none")

list1[[2]]

# plot 2
list2 <- prep_460_allsp_relabund(samples, reads, totals, 2000, 460, 0.01)
df <- list2[[1]]
(plot2 <- ggplot(data = df[df$df == 'NMDS',], aes(x = MDS1, y = MDS2)) + 
  geom_point(aes(color = group), size = 0.5, alpha = 0.5) +
  geom_path(data = df[df$df == "df_ell",], aes(x = MDS1, y = MDS2, colour = group), 
            size = 0.5) +
  geom_text(data = df[df$df == "NMDS.mean" & df$group == "K-Rat",], 
            aes(x = MDS1, y = MDS2,  
                label = .data$group[1], 
                color = .data$group[1]),
            size = 2) +
  geom_text(data = df[df$df == "NMDS.mean" & df$group == "CP: Control",], 
            aes(x = MDS1, y = MDS2,  
                label = .data$group[1], 
                color = .data$group[1]),
            size = 2) +
  geom_text(data = df[df$df == "NMDS.mean" & df$group == "CP: KR Exclosure",],
            aes(x = MDS1, y = MDS2,  
                label = .data$group[1], 
                color = .data$group[1]),
            size = 2) +
  geom_text(data = df %>% select(F.model, pval, min_total, min_rel_abund) %>% distinct(), 
            aes(x = Inf, y = Inf,
                label = paste("F.model = ", round(.data$F.model, 2),
                              "\n p = ", round(.data$pval, 4))),
            hjust = 1.1, vjust= 1.2, size = 2) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  labs(tag = "B") +  
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())) 

list2[[2]]

# plot 3
list3 <- prep_466_allsp_relabund(samples, reads, totals, 2000, 466, 0.01)
df <- list3[[1]]
(plot3 <- ggplot(data = df[df$df == 'NMDS',], aes(x = MDS1, y = MDS2)) + 
    geom_point(aes(color = group), size = 0.5, alpha = 0.5) +
    geom_path(data = df[df$df == "df_ell",], aes(x = MDS1, y = MDS2, colour = group), 
              size = 0.5) +
    # geom_text(data = df[df$df == "NMDS.mean" & df$group == "K-Rat",], 
    #           aes(x = MDS1, y = MDS2,  
    #               label = .data$group[1], 
    #               color = .data$group[1]),
    #           size = 2) +
    # geom_text(data = df[df$df == "NMDS.mean" & df$group == "CP: Control",], 
    #           aes(x = MDS1, y = MDS2,  
    #               label = .data$group[1], 
    #               color = .data$group[1]),
    #           size = 2) +
    # geom_text(data = df[df$df == "NMDS.mean" & df$group == "CP: KR Exclosure",],
    #           aes(x = MDS1, y = MDS2,  
    #               label = .data$group[1], 
    #               color = .data$group[1]),
    #           size = 2) +
    geom_text(data = list4[[2]], aes(x = -Inf, y = -Inf,
                                     label = paste0("psudeo F = ", round(list3[[2]]$F.Model[1], 2),
                                                    ", R^2 = ", round(list3[[2]]$R2[1], 2),
                                                    ", p = ", round(list3[[2]]$P.value[1], 3))),
              hjust = -0.05, vjust = -4, size = 3) +
    geom_text(data = list4[[2]], aes(x = -Inf, y = -Inf,
                                     label = paste0("psudeo F = ", round(list3[[2]]$F.Model[2], 2),
                                                    ", R^2 = ", round(list3[[2]]$R2[2], 2),
                                                    ", p = ", round(list3[[2]]$P.value[2], 3))),
              hjust = -0.05, vjust = -2.25, size = 3) +
    geom_text(data = list4[[2]], aes(x = -Inf, y = -Inf,
                                     label = paste0("KR \u2194 CP:C",
                                                    ", psudeo F = ", round(list3[[2]]$F.Model[3], 2),
                                                    ", R^2 = ", round(list3[[2]]$R2[3], 2),
                                                    ", p = ", round(list3[[2]]$P.value[3], 3))),
              hjust = -0.03, vjust = -0.75, size = 3) + 
    scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
    labs(tag = "C") +  
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank())) 

list3[[2]]$F.Model[2]

# plot 4
list4 <- prep_454_PPonly_relabund(samples_PP, reads, totals, 2000, 454, 0.01)
df <- list4[[1]]
(plot4 <- ggplot(data = df[df$df == 'NMDS',], aes(x = MDS1, y = MDS2)) + 
    geom_point(aes(color = group), size = 0.5, alpha = 0.5) +
    geom_path(data = df[df$df == "df_ell",], aes(x = MDS1, y = MDS2, colour = group), 
              size = 0.5) +
    # geom_text(data = df[df$df == "NMDS.mean" & df$group == "CP: Control",], 
    #           aes(x = MDS1, y = MDS2,  
    #               label = .data$group[1], 
    #               color = .data$group[1]),
    #           size = 2) +
    # geom_text(data = df[df$df == "NMDS.mean" & df$group == "CP: KR Exclosure",],
    #           aes(x = MDS1, y = MDS2,  
    #               label = .data$group[1], 
    #               color = .data$group[1]),
    #           size = 2) +
    geom_text(data = list4[[2]], aes(x = -Inf, y = -Inf,
                  label = paste0("psudeo F = ", round(list4[[2]]$F.Model, 2),
                                ", R^2 = ", round(list4[[2]]$R2, 2),
                                ", p = ", round(list4[[2]]$P.value, 3))),
              hjust = -0.05, vjust = -0.75, size = 3) +
    scale_color_manual(values = c("#E69F00", "#56B4E9")) +
    labs(tag = "D") +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())) 

round(list4[[2]]$F.Model, 2)

# plot 5
list5 <- prep_460_PPonly_relabund(samples_PP, reads, totals, 2000, 460, 0.01)
df <- list5[[1]]
(plot5 <- ggplot(data = df[df$df == 'NMDS',], aes(x = MDS1, y = MDS2)) + 
    geom_point(aes(color = group), size = 0.5, alpha = 0.5) +
    geom_path(data = df[df$df == "df_ell",], aes(x = MDS1, y = MDS2, colour = group), 
              size = 0.5) +
    # geom_text(data = df[df$df == "NMDS.mean" & df$group == "CP: Control",], 
    #           aes(x = MDS1, y = MDS2,  
    #               label = .data$group[1], 
    #               color = .data$group[1]),
    #           size = 2) +
    # geom_text(data = df[df$df == "NMDS.mean" & df$group == "CP: KR Exclosure",],
    #           aes(x = MDS1, y = MDS2,  
    #               label = .data$group[1], 
    #               color = .data$group[1]),
    #           size = 2) +
    geom_text(data = df %>% select(F.model, pval, min_total, min_rel_abund) %>% distinct(), 
              aes(x = -Inf, y = -Inf,
                  label = paste0("psudeo F = ", round(list5[[2]]$F.Model, 2),
                                 ", R^2 = ", round(list5[[2]]$R2, 2),
                                 ", p = ", round(list5[[2]]$P.value, 3))),
              hjust = -0.05, vjust = -0.75, size = 3) +
    scale_color_manual(values = c("#E69F00", "#56B4E9")) +
    labs(tag = "E") +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank())) 

list5[[2]]

# plot 6
list6 <- prep_466_PPonly_relabund(samples_PP, reads, totals, 2000, 466, 0.01)
df <- list6[[1]]
(plot6 <- ggplot(data = df[df$df == 'NMDS',], aes(x = MDS1, y = MDS2)) + 
    geom_point(aes(color = group), size = 0.5, alpha = 0.5) +
    geom_path(data = df[df$df == "df_ell",], aes(x = MDS1, y = MDS2, colour = group), 
              size = 0.5) +
    # geom_text(data = df[df$df == "NMDS.mean" & df$group == "CP: Control",], 
    #           aes(x = MDS1, y = MDS2,  
    #               label = .data$group[1], 
    #               color = .data$group[1]),
    #           size = 2) +
    # geom_text(data = df[df$df == "NMDS.mean" & df$group == "CP: KR Exclosure",],
    #           aes(x = MDS1, y = MDS2,  
    #               label = .data$group[1], 
    #               color = .data$group[1]),
    #           size = 2) +
    geom_text(data = df %>% select(F.model, pval, min_total, min_rel_abund) %>% distinct(), 
              aes(x = -Inf, y = -Inf,
                  label = paste0("psudeo F = ", round(list5[[2]]$F.Model, 2),
                                 ", R^2 = ", round(list5[[2]]$R2, 2),
                                 ", p = ", round(list5[[2]]$P.value, 3))),
              hjust = -0.05, vjust = -0.75, size = 3) +
    scale_color_manual(values = c("#E69F00", "#56B4E9")) +
    labs(tag = "F") +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank())) 

list6[[2]]

# multipanel plot

col1 <- ggplot() + annotate(geom = 'text', x = 1, y = 1, 
                            label = "Fall 2016", size = 5, fontface = "bold") + 
  theme_void() 
col2 <- ggplot() + annotate(geom = 'text', x = 1, y = 1, 
                            label = "Spring 2017", size = 5, fontface = "bold") + 
  theme_void() 
col3 <- ggplot() + annotate(geom = 'text', x = 1, y = 1, 
                            label = "Fall 2017", size = 5, fontface = "bold") + 
  theme_void() 


layoutplot <- "
cccdddiii##
eeefffaaa##
eeefffaaajj
ggghhhbbbjj
ggghhhbbb##
"
plot1 <- plot1 + theme(plot.margin = unit(c(0,20,0,0), "pt"))
plot2 <- plot2 + theme(plot.margin = unit(c(0,20,0,0), "pt"))

plotlist1 <- list(c = col1, d = col2, i = col3,
                 e = plot1, f = plot2, a = plot3,
                 g = plot4, h = plot5, b = plot6, j = legend)
wrap_plots(plotlist1, design = layoutplot) 

#ggsave("Plots/trnL_OTUs.png")

list1[[2]]
list2[[2]]
list3[[2]]
list4[[2]]
list5[[2]]
list6[[2]]

# ITS2 OTUS #-------------------------------------------------------------------

reads <- read_csv("Data/SequencedData/Plants/ProcessedData/ITS2_reads_WeeTU.csv")
totals <- read_csv("Data/SequencedData/Plants/ProcessedData/ITS2_totals.csv")
samples <- read_csv("Data/CollectionData/fecal_sample_collection.csv")

# get plant OTUs only
reads <- filter(reads, WTU.clade1 == 4) %>% 
  select(OTU:DataFrame)

# only PPs
samples_PP <- samples %>% filter(species == 'PP')

# plot 7
list7 <- prep_454_allsp_relabund_ITS2(samples, reads, totals, 2000, 454, 0.01)
df <- list7[[1]]
(plot7 <- ggplot(data = df[df$df == 'NMDS',], aes(x = MDS1, y = MDS2)) + 
    geom_point(aes(color = group), size = 0.5, alpha = 0.5) +
    geom_path(data = df[df$df == "df_ell",], aes(x = MDS1, y = MDS2, colour = group), 
              size = 0.5) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "K-Rat",], 
              aes(x = MDS1, y = MDS2,  
                  label = .data$group[1], 
                  color = .data$group[1]),
              size = 2) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "CP: Control",], 
              aes(x = MDS1, y = MDS2,  
                  label = .data$group[1], 
                  color = .data$group[1]),
              size = 2) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "CP: KR Exclosure",],
              aes(x = MDS1, y = MDS2,  
                  label = .data$group[1], 
                  color = .data$group[1]),
              size = 2) +
    geom_text(data = df %>% select(F.model, pval, min_total, min_rel_abund) %>% distinct(), 
              aes(x = Inf, y = Inf,
                  label = paste("F.model = ", round(.data$F.model, 2),
                                "\n p = ", round(.data$pval, 4))),
              hjust = 1.1, vjust= 1.2, size = 2) +
    scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
    labs(tag = "A") +
    theme_bw() +
    theme(legend.title = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank())) 

legend <- cowplot::get_legend(plot7)
plot7 <- plot7 + theme(legend.position = "none")

list7[[2]]

# plot 8
list8 <- prep_460_allsp_relabund_ITS2(samples, reads, totals, 2000, 460, 0.01)
df <- list8[[1]]
(plot8 <- ggplot(data = df[df$df == 'NMDS',], aes(x = MDS1, y = MDS2)) + 
    geom_point(aes(color = group), size = 0.5, alpha = 0.5) +
    geom_path(data = df[df$df == "df_ell",], aes(x = MDS1, y = MDS2, colour = group), 
              size = 0.5) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "K-Rat",], 
              aes(x = MDS1, y = MDS2,  
                  label = .data$group[1], 
                  color = .data$group[1]),
              size = 2) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "CP: Control",], 
              aes(x = MDS1, y = MDS2,  
                  label = .data$group[1], 
                  color = .data$group[1]),
              size = 2) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "CP: KR Exclosure",],
              aes(x = MDS1, y = MDS2,  
                  label = .data$group[1], 
                  color = .data$group[1]),
              size = 2) +
    geom_text(data = df %>% select(F.model, pval, min_total, min_rel_abund) %>% distinct(), 
              aes(x = Inf, y = Inf,
                  label = paste("F.model = ", round(.data$F.model, 2),
                                "\n p = ", round(.data$pval, 4))),
              hjust = 1.1, vjust= 1.2, size = 2) +
    scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
    labs(tag = "B") +  
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank())) 

list8[[2]]

# plot 9
list9 <- prep_466_allsp_relabund_ITS2(samples, reads, totals, 2000, 466, 0.01)
df <- list9[[1]]
(plot9 <- ggplot(data = df[df$df == 'NMDS',], aes(x = MDS1, y = MDS2)) + 
    geom_point(aes(color = group), size = 0.5, alpha = 0.5) +
    geom_path(data = df[df$df == "df_ell",], aes(x = MDS1, y = MDS2, colour = group), 
              size = 0.5) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "K-Rat",], 
              aes(x = MDS1, y = MDS2,  
                  label = .data$group[1], 
                  color = .data$group[1]),
              size = 2) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "CP: Control",], 
              aes(x = MDS1, y = MDS2,  
                  label = .data$group[1], 
                  color = .data$group[1]),
              size = 2) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "CP: KR Exclosure",],
              aes(x = MDS1, y = MDS2,  
                  label = .data$group[1], 
                  color = .data$group[1]),
              size = 2) +
    geom_text(data = df %>% select(F.model, pval, min_total, min_rel_abund) %>% distinct(), 
              aes(x = Inf, y = Inf,
                  label = paste("F.model = ", round(.data$F.model, 2),
                                "\n p = ", round(.data$pval, 4))),
              hjust = 1.1, vjust= 1.2, size = 2) +
    scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
    labs(tag = "C") +  
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank())) 

list9[[2]]

# plot 10
list10 <- prep_454_PPonly_relabund_ITS2(samples_PP, reads, totals, 2000, 454, 0.01)
df <- list10[[1]]
(plot10 <- ggplot(data = df[df$df == 'NMDS',], aes(x = MDS1, y = MDS2)) + 
    geom_point(aes(color = group), size = 0.5, alpha = 0.5) +
    geom_path(data = df[df$df == "df_ell",], aes(x = MDS1, y = MDS2, colour = group), 
              size = 0.5) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "CP: Control",], 
              aes(x = MDS1, y = MDS2,  
                  label = .data$group[1], 
                  color = .data$group[1]),
              size = 2) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "CP: KR Exclosure",],
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
          panel.grid.minor = element_blank())) 

list10[[2]]

# plot 11
list11 <- prep_460_PPonly_relabund_ITS2(samples_PP, reads, totals, 2000, 460, 0.01)
df <- list11[[1]]
(plot11 <- ggplot(data = df[df$df == 'NMDS',], aes(x = MDS1, y = MDS2)) + 
    geom_point(aes(color = group), size = 0.5, alpha = 0.5) +
    geom_path(data = df[df$df == "df_ell",], aes(x = MDS1, y = MDS2, colour = group), 
              size = 0.5) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "CP: Control",], 
              aes(x = MDS1, y = MDS2,  
                  label = .data$group[1], 
                  color = .data$group[1]),
              size = 2) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "CP: KR Exclosure",],
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
    labs(tag = "E") +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank())) 

list11[[2]]

# plot 6
list12 <- prep_466_PPonly_relabund_ITS2(samples_PP, reads, totals, 2000, 466, 0.01)
df <- list12[[1]]
(plot12 <- ggplot(data = df[df$df == 'NMDS',], aes(x = MDS1, y = MDS2)) + 
    geom_point(aes(color = group), size = 0.5, alpha = 0.5) +
    geom_path(data = df[df$df == "df_ell",], aes(x = MDS1, y = MDS2, colour = group), 
              size = 0.5) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "CP: Control",], 
              aes(x = MDS1, y = MDS2,  
                  label = .data$group[1], 
                  color = .data$group[1]),
              size = 2) +
    geom_text(data = df[df$df == "NMDS.mean" & df$group == "CP: KR Exclosure",],
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
    labs(tag = "F") +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank())) 

list12[[2]]

# multipanel plot

col1 <- ggplot() + annotate(geom = 'text', x = 1, y = 1, 
                            label = "Fall 2016", size = 5, fontface = "bold") + 
  theme_void() 
col2 <- ggplot() + annotate(geom = 'text', x = 1, y = 1, 
                            label = "Spring 2017", size = 5, fontface = "bold") + 
  theme_void() 
col3 <- ggplot() + annotate(geom = 'text', x = 1, y = 1, 
                            label = "Fall 2017", size = 5, fontface = "bold") + 
  theme_void() 


layoutplot <- "
cccdddiii##
eeefffaaa##
eeefffaaajj
ggghhhbbbjj
ggghhhbbb##
"
plot7 <- plot7 + theme(plot.margin = unit(c(0,20,0,0), "pt"))
plot8 <- plot8 + theme(plot.margin = unit(c(0,20,0,0), "pt"))

plotlist2 <- list(c = col1, d = col2, i = col3,
                 e = plot7, f = plot8, a = plot9,
                 g = plot10, h = plot11, b = plot12, j = legend)
wrap_plots(plotlist2, design = layoutplot) 

#ggsave("Plots/ITS2_OTUs.png")

list7[[2]]
list8[[2]]
list9[[2]]
list10[[2]]
list11[[2]]
list12[[2]]
