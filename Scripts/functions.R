# FUNCTIONS #
# 19 Nov 2019
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

# convert trnL from reads to proportions #
