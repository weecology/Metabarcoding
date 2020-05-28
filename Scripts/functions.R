# FUNCTIONS #
# May 2020
# EKB

source('Scripts/find_millet_OTUs.R')

read_in_trnL_files <- function(path = "Data/SequencedData/Plants/RawData/PreppedFiles/trnL",
                               files = dir(path, pattern = "*.csv")){
  # this function reads all CSV files from set path and puts them in a list
  data <- files %>% 
    map(~ read_csv(file.path(path, .)))
  return(data)
}


read_in_ITS2_files <- function(path = "Data/SequencedData/Plants/RawData/PreppedFiles/ITS2",
                               files = dir(path, pattern = "*.csv")){
  # this function reads all CSV files from set path and puts them in a list
  data <- files %>% 
    map(~ read_csv(file.path(path, .)))
  return(data)
}


add_plot_type <- function(data){
  # add plot type to dataset -- control or KR exclosure
  data$plot_type <- NA
  for (i in 1:length(data$plot_type)) {
    if (data$plot[i] %in% c(4, 11, 14, 17)) {
      data$plot_type[i] = 'Control'
    } else {
      data$plot_type[i] = 'KR_Exclosure'
    }
  }
  return(data)
}


add_plotting_group <- function(data){
  # create grouping column based on species and plot
  data$group = NA                     
  for (i in 1:length(data$species)) {
    if (data$species[i] %in% c('DO', 'DM')) {
      data$group[i] = 'Krat'
    } else if (data$plot[i] %in% c(4, 11, 14, 17)) {
      data$group[i] = 'PP_control'
    } else {
      data$group[i] = 'PP_exclosure'
    }
  }
  return(data)
}


summarize_trnL_by_WeeTU <- function(data, col_quotes, col_no_quotes){
  
  # make dataframe with all unique species
  all_species <- data %>% 
    subset(., !duplicated(WTU.species)) %>% 
    select(Kingdom:WTU.species)
  
  # select unique values to specified taxa level
  if (col_quotes %in% c("Species", "WTU.species")) {
    all_taxa <- all_species
  } else if (col_quotes %in% c("Genus", "WTU.genus")) {
    all_taxa <- all_species %>% select(-Species, -WTU.species)
  } else if (col_quotes %in% c("Subfamily", "WTU.subfamily")) {
    all_taxa <- all_species %>% select(-Species, -WTU.species,
                                       -Genus, -WTU.genus)
  } else if (col_quotes %in% c("Family", "WTU.family")) {
    all_taxa <- all_species %>% select(-Species, -WTU.species,
                                       -Genus, -WTU.genus,
                                       -Subfamily, -WTU.subfamily)
  } else if (col_quotes %in% c("Order", "WTU.order")) {
    all_taxa <- all_species %>% select(-Species, -WTU.species,
                                       -Genus, -WTU.genus,
                                       -Subfamily, WTU.subfamily,
                                       -Family, -WTU.family)
  } else if (col_quotes %in% c("Clade2", "WTU.clade2")) {
    all_taxa <- all_species %>% select(Kingdom, WTU.kingdom,
                                       Clade1, WTU.clade1,
                                       Clade2, WTU.clade2)
  } else if (col_quotes %in% c("Clade1", "WTU.clade1")) {
    all_taxa <- all_species %>% select(Kingdom, WTU.kingdom, Clade1, WTU.clade1)
  } else  {
    all_taxa <- all_species %>% select(Kingdom, WTU.kingdom)
  }
  
  column <- enquo(col_no_quotes)
  
  sum <- data %>%
    group_by(SampleID, !! column) %>%
    summarise(Reads = sum(Reads))
  sum <- left_join(sum, all_taxa) %>% distinct()
  
  return(sum)
  
}


summarize_ITS2_by_WeeTU <- function(data, col_quotes, col_no_quotes){
  
  # make dataframe with all unique species
  all_species <- data %>% 
    subset(., !duplicated(WTU.species)) %>% 
    select(Domain:WTU.species)
  
  # select unique values to specified taxa level
  if (col_quotes %in% c("Species", "WTU.species")) {
    all_taxa <- all_species
  } else if (col_quotes %in% c("Genus", "WTU.genus")) {
    all_taxa <- all_species %>% select(-Species, -WTU.species)
  } else if (col_quotes %in% c("Family", "WTU.family")) {
    all_taxa <- all_species %>% select(-Species, -WTU.species,
                                       -Genus, -WTU.genus)
  } else if (col_quotes %in% c("Order", "WTU.order")) {
    all_taxa <- all_species %>% select(-Species, -WTU.species,
                                       -Genus, -WTU.genus,
                                       - Family, -WTU.family)
  } else if (col_quotes %in% c("Class", "WTU.class")) {
    all_taxa <- all_species %>% select(-Species, -WTU.species,
                                       -Genus, -WTU.genus,
                                       -Family, -WTU.family,
                                       -Order, -WTU.order)
  } else if (col_quotes %in% c("Clade1", "WTU.clade1")) {
    all_taxa <- all_species %>% select(Domain, WTU.domain, Clade1, WTU.clade1)
  } else  {
    all_taxa <- all_species %>% select(Domain, WTU.domain)
  }
  
  column <- enquo(col_no_quotes)
  
  sum <- data %>%
    group_by(SampleID, !! column) %>%
    summarise(Reads = sum(Reads))
  sum <- left_join(sum, all_taxa) %>% distinct()
  
  return(sum)
  
}


filter_reads_data_trnL <- function(samples,
                              reads,
                              totals,
                              yr = c(2016, 2017),
                              reads_min = 2000,
                              rel_reads_min = 0.001){
  
  # add plot type to fecal collection data
  # add group for plotting
  # and remove samples that were part of the trap/bait test
  samples <- add_plot_type(samples) %>% 
    add_plotting_group() %>% 
    filter(is.na(notes), year %in% yr) 
  
  # select only fecal samples
  fecal_id <- samples$vial_barcode
  
  # add totals to reads df
  # select only fecal samples and remove millet OTUs
  # and make relative reads column
  reads <- full_join(reads, totals)
  reads <- reads %>% 
    filter(SampleID %in% fecal_id, !OTU %in% millet_OTUs$OTU) %>% 
    mutate(Rel_Reads = Reads/Total_Reads)
  
  # filter data by minimum total reads and/or minimum relative reads
  reads <- reads %>% 
    filter(Total_Reads >= reads_min, Rel_Reads >= rel_reads_min)
  
  return_list <- list(samples, fecal_id, reads)
  names(return_list) <- c("samples", "fecal_id", "reads")
  return(return_list)
  
}

filter_reads_data_ITS2 <- function(samples,
                                   reads,
                                   totals,
                                   yr = c(2016, 2017),
                                   reads_min = 2000,
                                   rel_reads_min = 0.001){
  
  # add plot type to fecal collection data
  # add group for plotting
  # and remove samples that were part of the trap/bait test
  samples <- add_plot_type(samples) %>% 
    add_plotting_group() %>% 
    filter(is.na(notes), year %in% yr) 
  
  # select only fecal samples
  fecal_id <- samples$vial_barcode
  
  # add totals to reads df
  # select only fecal samples and remove millet OTUs
  # and make relative reads column
  reads <- full_join(reads, totals)
  reads <- reads %>% 
    filter(SampleID %in% fecal_id, !OTU %in% millet_OTUs_ITS2_no.hirt$OTU) %>% 
    mutate(Rel_Reads = Reads/Total_Reads)
  
  # filter data by minimum total reads and/or minimum relative reads
  reads <- reads %>% 
    filter(Total_Reads >= reads_min, Rel_Reads >= rel_reads_min)
  
  return_list <- list(samples, fecal_id, reads)
  names(return_list) <- c("samples", "fecal_id", "reads")
  return(return_list)
  
}


data_prep_multivariate <- function(list){
  
  # this function preps data for running NMDS
  # list argument is output from `filter_reads_data` fxn
  
  samples <- list[[1]]
  reads <- list[[3]]
  
  # convert to appropriate format for NMDS
  reads <- select(reads, SampleID, OTU, Rel_Reads)
  reads_spread <- pivot_wider(reads, names_from = OTU, values_from = Rel_Reads)
  reads_spread[is.na(reads_spread)] = 0
  reads_spread <- reads_spread %>% tibble::column_to_rownames("SampleID")
  
  sampleID <- intersect(reads$SampleID, samples$vial_barcode)
  groups <- samples %>% 
    filter(vial_barcode %in% sampleID)
  groups <- groups[match(sampleID, groups$vial_barcode),]
  
  return_list <- list(reads_spread, sampleID, groups)
  names(return_list) <- c("reads_spread", "sampleID", "groups")
  return(return_list)
  
}


run_metaMDS_til_converge <- function(df, prev_best, dist_metric, n) {
  repeat {
    # do something
    prev_best <- metaMDS(df, distance = dist_metric, 
                         trymax = 50, k = n, 
                         previous.best = prev_best)
    # exit if condition is met
    if (prev_best$converged == TRUE) break
  }
  return(prev_best)
}


binarize <- function(df){
  df[df > 0] <- 1
  return(df)
}


NMDS_plotting_prep <- function(data_list, dist_matrix) {
  
  # this function makes ellipses for plotting and gets perMANOVA values
  # data_list = output of `data_prep_multivariate` fxn
  
  # make dataframe with with MDS values and grouping variable
  groups <- data_list[[3]]
  NMDS <- data.frame(MDS1 = dist_matrix$points[,1], 
                     MDS2 = dist_matrix$points[,2], 
                     group = groups$group)
  
  # get mean point for each group
  NMDS.mean <- aggregate(NMDS[,1:2], list(group = groups$group), mean)
  
  # save results of ordiellipse() as an object
  plot(dist_matrix$points)
  ord <- ordiellipse(dist_matrix, 
                     groups$group, 
                     display = "sites", 
                     kind = "se", 
                     conf = 0.95, 
                     label = T)
  
  # make ellipses dataframe for ggplot2 
  df_ell <- data.frame()
  for(g in levels(NMDS$group)) {
    df_ell <-
      rbind(df_ell, cbind(as.data.frame(with(
        NMDS[NMDS$group == g,],
        vegan:::veganCovEllipse(ord[[g]]$cov, ord[[g]]$center, ord[[g]]$scale)
      )),
      group = g))
  }
  
  # run perMANOVA
  group = as.matrix(groups$group)
  perMANOVA_output <- adonis(data_list[[1]] ~ group, permutations = 10000)
  
  # list of objects to return
  return_list <- list(NMDS, NMDS.mean, df_ell, perMANOVA_output)
  names(return_list) <- c("NMDS", "NMDS.mean", "df_ell", "perMANOVA_output")
  return(return_list)
  
}


plot_NMDS_ggplot2 <- function(NMDS_list) {
  
  # NMDS_list is output list from `NMDS_plotting_prep`
  
  # ggplot of NMDS
  plot <- ggplot(data = NMDS_list[[1]], aes(x = MDS1, y = MDS2)) + 
    geom_point(aes(color = group)) +
    geom_path(data = NMDS_list[[3]], aes(x = NMDS1, y = NMDS2, colour = group), 
              size = 1) +
    geom_text(aes(x = NMDS_list[[2]]$MDS1[1], y = NMDS_list[[2]]$MDS2[1], 
                  label = NMDS_list[[2]]$group[1], color = NMDS_list[[2]]$group[1])) +
    geom_text(aes(x = NMDS_list[[2]]$MDS1[2], y = NMDS_list[[2]]$MDS2[2], 
                  label = NMDS_list[[2]]$group[2], color = NMDS_list[[2]]$group[2])) +
    geom_text(aes(x = NMDS_list[[2]]$MDS1[3], y = NMDS_list[[2]]$MDS2[3], 
                  label = NMDS_list[[2]]$group[3], color = NMDS_list[[2]]$group[3])) +
    scale_color_manual(values = cbPalette) +
    theme_bw() +
    theme(legend.position = 'none',
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    annotate(geom = "text", x = Inf, y = Inf, hjust = 1.1, vjust= 1.2,
             label = paste("atop(' F.model = '*", round(NMDS_list[[4]]$aov.tab$F.Model[1], 2),"
                         ,' p = '*", round(NMDS_list[[4]]$aov.tab$`Pr(>F)`[1], 4),")"), parse=T)
  
  return(plot)
  
}
