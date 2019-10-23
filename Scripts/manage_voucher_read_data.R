library(dplyr)
# Takes data from ITS and TrnL voucher data files, calculates total
# number of reads for each sample and adds that to the voucher collection file
# Data comes from repo.

add_all_voucher_readtotals = function(){
### Checking current state of trnL collection files

trnL_data = read.csv("./SequencedData/Plants/trnL_voucher_data.csv")
ITS_data = read.csv("./SequencedData/Plants/weeTU_ITS_voucher_data.csv")
collection = read.csv("./CollectionData/plant_voucher_collection.csv")

# Add trnL total reads to collection file

sample_trnL_reads = trnL_data %>% group_by(Sample) %>% summarise(trnl_total=sum(Reads))
sample_ITS_reads = ITS_data %>% group_by(Sample) %>% summarise(ITS_total=sum(Reads))
add_trnL_total = left_join(collection,sample_trnL_reads, by=c("vial_barcode" = "Sample"))
add_ITS_total = left_join(add_trnL_total, sample_ITS_reads, by=c("vial_barcode" = "Sample"))
write.csv(add_ITS_total, "./CollectionData/plant_voucher_collection.csv")
}

