# Clean Fecal Collection Data from 460 and (-)460
# EKB
# 4/25/2017

#=================================================
# LIBRARIES

library(dplyr)
library(sqldf)

#=================================================
# FUNCTIONS

compare_tags = function(ws,scannerfile) {
  
  # extract tag numbers from raw data
  sheets = subset(ws$tag,!is.na(ws$tag))
  
  # load data from scanner
  scandat = read.table(scannerfile, 
                       header=FALSE, 
                       sep='.', 
                       blank.lines.skip=TRUE,
                       col.names=c('v1','tag','date','time'))
  
  # extract 6-digit tag numbers from scanner
  scans = vector()
  for (tag in as.vector(scandat$tag)) {
    scans = append(scans,substr(tag,5,10))
  }
  
  scannotsheet = setdiff(scans,sheets)
  sheetnotscan = setdiff(sheets,scans)
  unpaired = data.frame(where=c(rep('scan',length(scannotsheet)),rep('sheet',length(sheetnotscan))),
                        tag=c(scannotsheet,sheetnotscan))
  
  print(unpaired)
}

#=================================================
# PERIOD (-)460 #

### get files ###

newperiod = '-460'
filepath = 'C:/Users/ellen.bledsoe/Dropbox/Portal/PORTAL_primary_data/DNA/Rodents/'
scannerfile = paste(filepath, 'tags', newperiod, '.txt', sep='')

ws <- read.csv('C:/Users/ellen.bledsoe/Dropbox/Portal/PORTAL_primary_data/DNA/Rodents/fecal_samples.csv', header = TRUE)
ws <- ws %>% filter(period == '-460') %>% rename(tag = PIT_tag)

### check tags ###

compare_tags(ws, scannerfile)
#     - only problem are tags that are on sheet and not scans, unless sp code with number

### check characteristics ###

olddat = read.csv('C:/Users/ellen.bledsoe/Desktop/Git/PortalData/Rodents/Portal_rodent.csv',na.strings='',as.is=T)

# Subset of most recent four years of data, for comparing recaptures
recentdat = olddat[olddat$yr >= as.numeric(ws$yr[1])-3,]

# Check sex/species on recaptures
#    -conflicts can be resolved if there's a clear majority, or if clear sexual characteristics
#    -also look back in book to see if sex/species data was manually changed before for a particular tag number
#    -when making changes to old or new data, note in book
sqldf("SELECT recentdat.period, recentdat.plot, ws.plot, recentdat.species, ws.species, recentdat.sex, ws.sex, recentdat.tag
      FROM recentdat INNER JOIN ws ON recentdat.tag = ws.tag
      WHERE (((recentdat.species)<>(ws.species)) And ((recentdat.tag)=(ws.tag))) Or (((recentdat.sex)<>(ws.sex)));")

#===================================================
# PERIOD 460 #

### get files ###

newperiod = '460'
filepath = 'C:/Users/ellen.bledsoe/Dropbox/Portal/PORTAL_primary_data/Rodent/Raw_data/New_data/'
scannerfile = paste(filepath, 'tag scans/tags', newperiod, '.txt', sep='')

ws <- read.csv('C:/Users/ellen.bledsoe/Dropbox/Portal/PORTAL_primary_data/DNA/Rodents/fecal_samples.csv', header = TRUE)
ws <- ws %>% filter(period == '460') %>% rename(tag = PIT_tag)

### check tags ###

compare_tags(ws, scannerfile)
#     - only problem are tags that are on sheet and not scans, unless sp code with number

### check characteristics ###

olddat = read.csv('C:/Users/ellen.bledsoe/Desktop/Git/PortalData/Rodents/Portal_rodent.csv',na.strings='',as.is=T)

# Subset of most recent four years of data, for comparing recaptures
recentdat = olddat[olddat$yr >= as.numeric(ws$yr[1])-3,]

# Check sex/species on recaptures
#    -conflicts can be resolved if there's a clear majority, or if clear sexual characteristics
#    -also look back in book to see if sex/species data was manually changed before for a particular tag number
#    -when making changes to old or new data, note in book
sqldf("SELECT recentdat.period, recentdat.plot, ws.plot, recentdat.species, ws.species, recentdat.sex, ws.sex, recentdat.tag
      FROM recentdat INNER JOIN ws ON recentdat.tag = ws.tag
      WHERE (((recentdat.species)<>(ws.species)) And ((recentdat.tag)=(ws.tag))) Or (((recentdat.sex)<>(ws.sex)));")
