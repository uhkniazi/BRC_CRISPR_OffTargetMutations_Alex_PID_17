# File: 01_createDbEntry.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: list the samples and create appropriate db entries
# Date: 05/07/2019


## set variables and source libraries
source('header.R')

## connect to mysql database 
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

# sample and file table
dbGetQuery(db, paste('describe Sample;'))
cSampleCol = dbGetQuery(db, paste('describe Sample;'))$Field[-1]

dbGetQuery(db, paste('describe File;'))
cFileCol = dbGetQuery(db, paste('describe File;'))$Field[-1]

# setwd(gcRemoteDir)
setwd('dataExternal/')
setwd('raw/')
setwd('raw/S198/processed/')

# list the files
cvFiles = list.files(pattern = 'fastq.gz')
i = grepl('Uni', cvFiles) & grepl('R1_|R2_', cvFiles)
cvFiles = cvFiles[i]
# each sample has 2 files 
fSplit = gsub('^(Uni\\d)_.+', '\\1', cvFiles)

lFiles = split(cvFiles, fSplit)

setwd(gcswd)

## create the meta data table for sample table
cSampleCol

dfSamples = data.frame(idProject=g_pid, idData=g_did, title=unique(fSplit), 
                       description='S198 run for Alex CRISPR Cas9 sequencing for off target analysis',
                       group1 = unique(fSplit))

# write this data to the database
rownames(dfSamples) = NULL

### NOTE: Do not execute this anymore as entry created
# # write this table to database
# dbWriteTable(db, name='Sample', value=dfSamples, append=T, row.names=F)
# get this table again from database with ids added
g_did
dfSamples = dbGetQuery(db, paste('select * from Sample where Sample.idData = 35;'))

# create entries for these files in the database
dbListTables(db)
cn = dbListFields(db, 'File')[-1]
cn
# get the names of the samples
temp = lapply(dfSamples$title, function(x){
  # get the file names
  df = data.frame(name=lFiles[[x]], type='fastq', idSample=dfSamples[dfSamples$title == x, 'id'])
  return(df)
})

dfFiles = do.call(rbind, temp)
rownames(dfFiles) = NULL

# write this table to database
## note: do not execute as it is already done
# dbWriteTable(db, name='File', value=dfFiles, append=T, row.names=F)

dbDisconnect(db)