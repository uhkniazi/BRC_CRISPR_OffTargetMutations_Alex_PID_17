# File: 07_bam_files_qa.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: quality checks on the bam files
# Date: 9/10/2018


## set variables and source libraries
source('header.R')
library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CBamQuality/master/CBamQuality.R'
download(url, 'CBamQuality.R')

# load the required packages
source('CBamQuality.R')
# delete the file after source
unlink('CBamQuality.R')

## choose the sequence names over which to calculate coverage
cvSeqnames = c(paste0('chr', 1:22), 'chrX', 'chrY')

## connect to mysql database to get sample information
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'File')
dbListFields(db, 'Sample')
# get the query
g_did
q = paste0('select Sample.id as sid, Sample.group1, Sample.group2, Sample.group3, Sample.title, Sample.description, File.* from Sample, File
           where (Sample.idData = 35) AND (File.idSample = Sample.id AND File.type like "%sort%")')
dfSample = dbGetQuery(db, q)
nrow(dfSample)
dfSample
# close connection after getting data
dbDisconnect(db)
# remove any whitespace from the names
dfSample$name = gsub(" ", "", dfSample$name, fixed = T)
dfSample$title = gsub(" ", "", dfSample$title, fixed = T)

#### set working directory to appropriate location with bam files
setwd('dataExternal/raw/Aligned/')
csFiles = list.files('.', pattern = '*.bam$')
# check if these files match the file names in database
table(dfSample$name %in% csFiles)

dfSample = dfSample[dfSample$name %in% csFiles, ]
dim(dfSample)
csFiles = dfSample$name

## perform the analysis one sample at a time
## function to write the qa files
write.qa = function(bfn, sn){
  # open bam file
  bf = BamFile(bfn)
  lBam = lapply(sn, function(x){
    return(CBamScaffold(bf, x))
  })
  cat(paste('done', bfn, '\n'))
  return(lBam)
}

lAllBams = lapply(csFiles, function(x){
  write.qa(x, cvSeqnames)
})
names(lAllBams) = dfSample$id

setwd(gcswd)
n = make.names(paste('CBamScaffold S198 run for Alex CRISPR Cas9 sequencing for off target analysis data id 35 rds'))
lAllBams$meta = dfSample
lAllBams$desc = paste('CBamScaffold S198 run for Alex CRISPR Cas9 sequencing for off target analysis data id 35', date())
n2 = paste0('~/Data/MetaData/', n)
save(lAllBams, file=n2)

lAllBams

# comment out as this has been done once
# library('RMySQL')
# db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
# dbListTables(db)
# dbListFields(db, 'MetaFile')
# df = data.frame(idData=g_did, name=n, type='rds', location='~/Data/MetaData/',
#                 comment='CBamScaffold S198 run for Alex CRISPR Cas9 sequencing for off target analysis data id 35')
# dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
# dbDisconnect(db)

### create the plots of interest
getwd()
setwd(gcswd)
dfSample = lAllBams$meta
lAllBams$desc = NULL
lAllBams$meta = NULL
names(lAllBams) = dfSample$title


pdf(file='results/bam.qa.pdf')
par(mfrow=c(2,2))
## returns the coverage for each scaffold in the bam file in matrix
f1 = function(ob){
  mat = sapply(ob, CBamScaffold.getCoverage)
  n = sapply(ob, CBamScaffold.getSeqname)
  colnames(mat) = n
  return(log(mat+1))
}

## call the function f1
lMat = lapply(lAllBams, f1)
i = seq_along(cvSeqnames)

## reorder the matrix, so that same chromosome/scaffold from all bam files are in one matrix
lMat.ordered = lapply(i, function(x){
  m = sapply(lMat, function(y) y[,x])
})

names(lMat.ordered) = cvSeqnames
iCol = rainbow(length(lAllBams))
## plot the lowess profile 
temp = sapply(cvSeqnames, function(x){
  m = apply(lMat.ordered[[x]], 2, function(y) lowess(y, f=2/10)$y )
  matplot(m, type='l', main=paste('Lowess fit to coverage', x), col=iCol, xlab='Bins', sub='Binned Coverage', ylab='Coverage', 
          lty=1:length(iCol))
})

par(mfrow=c(1,1))
plot.new()
legend('center', legend = names(lAllBams), ncol=3, col = iCol, lty=1:length(iCol), cex=0.6, lwd=2)


## average coverage, read width
## modify function with replace = T 
CBamScaffold.getReadWidthSample = function(obj, size=1000) sample(obj@ivReadWidth, size = size, replace = T)

f1 = function(ob){
  mat = sapply(ob, CBamScaffold.getReadWidthSample)
  n = sapply(ob, CBamScaffold.getSeqname)
  colnames(mat) = n
  return(mean(colMeans(mat)))
}

## call the function f1
ivMat = sapply(lAllBams, f1)
barplot(ivMat, las=2, main='Average Read Width', ylab = 'Average Read Width', cex.names =0.6)

# number of reads aligned
f1 = function(ob){
  n = sapply(ob, function(x) length(CBamScaffold.getReadWidth(x)))
  n = sum(n)/1e+6
  return(n)
}

iReadCount = sapply(lAllBams, f1)

barplot(iReadCount, las=2, main='No. of reads aligned', ylab = 'No. of Reads in Millions', cex.names =0.8)

# average binned coverage distribution
f1 = function(ob){
  mat = sapply(ob, getCoverageGammaParam)['shape',]
  n = sapply(ob, CBamScaffold.getSeqname)
  names(mat) = n
  return(mat)
}

## call the function f1
lMat = lapply(lAllBams, function(x) tryCatch(expr = f1(x), error=function(e)NULL))

boxplot(lMat, las=2, main='Average Binned Coverage', ylab = 'Average', outline=F, xaxt='n')
axis(1, at=1:length(lMat), labels = names(lMat), cex.axis=0.7, las=2)

dev.off(dev.cur())

