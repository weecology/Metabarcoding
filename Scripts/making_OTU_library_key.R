# Making OTU Library Key
# 12 November 2019
# EKB

# This script is for making an OTU library for our plant vouchers. 
#   - there will be separate files for trnL and ITS2
#   - each file will contain:
#       - the OTU from Jonah Ventures
#       - the family, genus, and species of the OTU
#       - the WeeTU for filtering by family, genus, or species

# Current questions:
#   - Are we basing this only on plant vouchers that were sent in?
#   - OR from Portal plant list?
#   - OR are we going through OTUs from Jonah Ventures that make family, genus,
#     species that we have at the site?
#   - OR the OTUs that have number/proportion of reads above a certain threshold?

# TBD before really digging in: do I have data from all the plant vouchers?
#   - pull out voucher data first? Then convert to reads or proportions or not
#     even bother with that part?
#   - might want to source some code for pulling out sample numbers, so you can 
#     use new vial_id csv for pulling out plant voucher samples

#==============================================================================

# LIBRARIES & SOURCES #
library(tidyverse)
source("Scripts/make_trnL_SEED_dataframe.R")
