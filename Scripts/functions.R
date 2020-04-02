# FUNCTIONS #
# April 2020
# EKB

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


filter_reads_data <- function(samples,
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
  # select only fecal samples
  # and make relative reads column
  reads <- full_join(reads, totals)
  reads <- reads %>% 
    filter(SampleID %in% fecal_id) %>% 
    mutate(Rel_Reads = Reads/Total_Reads)
  
  # filter data by minimum total reads and/or minimum relative reads
  reads <- reads %>% 
    filter(Total_Reads >= reads_min, Rel_Reads >= rel_reads_min)
  
  return_list <- list(samples, fecal_id, reads)
  names(return_list) <- c("samples", "fecal_id", "reads")
  return(return_list)
  
}

data_prep_for_NMDS <- function(list){
  
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
  
  return_list <- list(reads_spread, sampleID, groups)
  names(return_list) <- c("reads_spread", "sampleID", "groups")
  return(return_list)
  
}
