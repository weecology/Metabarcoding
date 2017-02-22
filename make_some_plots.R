# Make plots with both ITS and trnL data
# Feb 22, 2017

########################
# SOURCES

setwd("C:/Users/ellen.bledsoe/Desktop/Git/Metagenomics")
source("./trnL_data_prep.R")
source("./ITS2_data_prep.R")
source("./quick_plotting_functions.R")

########################
# PLOTS

### trnL

# rank abundance 
rank_abundance(samples, reads, sp = c('PP', 'DM', 'DO')) 
rank_abundance(samples, reads, sp = 'PP')
rank_abundance(samples, reads, sp = c('DM', 'DO'))

# sum by family
sum_by_family(taxa_trnL, samples, reads)
sum_by_family(taxa_trnL, samples, reads, sp = 'PP')
sum_by_family(taxa_trnL, samples, reads, sp = c('DM', 'DO'))
sum_by_family(taxa_trnL, samples, reads, sp = 'DM')
sum_by_family(taxa_trnL, samples, reads, sp = 'DO')

# by individual
trap_vs_fresh_indv(samples, reads)


### ITS

# rank abundance 
rank_abundance(samples, all_ITS, cut_off = 25) 
rank_abundance(samples, all_ITS, sp = 'PP', cut_off = 25)
rank_abundance(samples, all_ITS, sp = c('DM', 'DO'), cut_off = 25)

# sum by family
sum_by_family(taxa_its, samples, all_ITS)
sum_by_family(taxa_its, samples, all_ITS, sp = 'PP')
sum_by_family(taxa_its, samples, all_ITS, sp = c('DM', 'DO'))
sum_by_family(taxa_its, samples, reads, sp = 'DM')
sum_by_family(taxa_its, samples, reads, sp = 'DO')

# by individual
trap_vs_fresh_indv(samples, all_ITS, cut_off = 100)
