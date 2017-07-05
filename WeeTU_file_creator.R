### Format ITS files to use WeeTU
library(readr)


convert_OTU_to_WeeTU = function(data, OTU_col_number){
  OTU_number = readr::parse_number(data[,OTU_col_number])
  WeeTU = paste(OTU_number, ".WeeTU", sep="")
  data = data[-OTU_col_number]
  data = cbind(WeeTU,data)
  return(data)
}
path = "./SequencedData/Plants/"
filenames = c("ITS_BLAST_taxa_link_file.csv", "ITS_fecal_data.csv", 
              "ITS_sequences_from_fna.csv", "ITS_voucher_data.csv")
for(file in filenames){
  data = read.csv(paste(path, file, sep=""))
  data = convert_OTU_to_WeeTU(data,1)
  write.csv(data,paste(path,"weetu_", file, sep=""), row.names=FALSE)
}



 
